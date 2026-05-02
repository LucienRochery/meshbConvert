#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <unordered_map>

extern "C" {
#include <libmeshb7.h>
}

// ---- Vertex deduplication via spatial hashing ----

struct Vec3Key {
    uint64_t bits[3];
    bool operator==(const Vec3Key &other) const
    {
        return bits[0] == other.bits[0] && bits[1] == other.bits[1] && bits[2] == other.bits[2];
    }
};

struct Vec3KeyHash {
    size_t operator()(const Vec3Key &kk) const
    {
        size_t hh = 14695981039346656037ULL;
        for (int ii = 0; ii < 3; ii++) {
            hh ^= kk.bits[ii];
            hh *= 1099511628211ULL;
        }
        return hh;
    }
};

static uint64_t edgeKey(int va, int vb)
{
    uint32_t lo = (uint32_t)(va < vb ? va : vb);
    uint32_t hi = (uint32_t)(va < vb ? vb : va);
    return ((uint64_t)hi << 32) | lo;
}

static Vec3Key makeKey(float xx, float yy, float zz)
{
    Vec3Key kk;
    kk.bits[0] = 0; kk.bits[1] = 0; kk.bits[2] = 0;
    memcpy(&kk.bits[0], &xx, sizeof(float));
    memcpy(&kk.bits[1], &yy, sizeof(float));
    memcpy(&kk.bits[2], &zz, sizeof(float));
    return kk;
}

static bool isBinaryStl(FILE *fp)
{
    char header[80];
    if (fread(header, 1, 80, fp) != 80) { rewind(fp); return false; }
    uint32_t nTri;
    if (fread(&nTri, 4, 1, fp) != 1) { rewind(fp); return false; }
    fseek(fp, 0, SEEK_END);
    long fileSize = ftell(fp);
    rewind(fp);
    long expectedBinary = 80 + 4 + (long)nTri * 50;
    return fileSize == expectedBinary;
}

static bool endsWith(const char *str, const char *suffix)
{
    int lenStr = (int)strlen(str);
    int lenSuf = (int)strlen(suffix);
    if (lenSuf > lenStr) return false;
    return strcmp(str + lenStr - lenSuf, suffix) == 0;
}

static void printUsage(const char *prog)
{
    printf("Usage: %s input output\n", prog);
    printf("\n");
    printf("Converts between STL and meshb surface mesh formats.\n");
    printf("Direction is detected from file extensions.\n");
    printf("\n");
    printf("  STL  -> meshb : %s model.stl model.meshb\n", prog);
    printf("  meshb -> STL  : %s model.meshb model.stl\n", prog);
    printf("\n");
    printf("STL to meshb:\n");
    printf("  - Reads binary or ASCII STL\n");
    printf("  - Deduplicates vertices\n");
    printf("  - Extracts boundary edges\n");
    printf("  - Assigns connected-component refs to triangles\n");
    printf("\n");
    printf("meshb to STL:\n");
    printf("  - Writes binary STL from meshb triangles\n");
    printf("  - Quads are split into two triangles\n");
    printf("  - Other element types (corners, edges, volumes) are dropped with a warning\n");
}

// ---- STL to meshb ----

static int stlToMeshb(const char *inputFile, const char *outputFile)
{
    FILE *fp = fopen(inputFile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    std::vector<double> vx, vy, vz;
    std::vector<int> t0, t1, t2;
    std::unordered_map<Vec3Key, int, Vec3KeyHash> vertexMap;

    auto addVertex = [&](float xx, float yy, float zz) -> int {
        Vec3Key kk = makeKey(xx, yy, zz);
        auto [it, inserted] = vertexMap.emplace(kk, (int)vx.size());
        if (inserted) {
            vx.push_back((double)xx);
            vy.push_back((double)yy);
            vz.push_back((double)zz);
        }
        return it->second;
    };

    if (isBinaryStl(fp)) {
        fseek(fp, 80, SEEK_SET);
        uint32_t nTri;
        fread(&nTri, 4, 1, fp);
        printf("Reading binary STL: %u triangles\n", nTri);
        for (uint32_t ii = 0; ii < nTri; ii++) {
            float buf[12];
            uint16_t attr;
            fread(buf, sizeof(float), 12, fp);
            fread(&attr, 2, 1, fp);
            int va = addVertex(buf[3], buf[4], buf[5]);
            int vb = addVertex(buf[6], buf[7], buf[8]);
            int vc = addVertex(buf[9], buf[10], buf[11]);
            t0.push_back(va);
            t1.push_back(vb);
            t2.push_back(vc);
        }
    } else {
        fseek(fp, 0, SEEK_SET);
        char line[256];
        float triVerts[9];
        int vCount = 0;
        while (fgets(line, sizeof(line), fp)) {
            float xx, yy, zz;
            if (sscanf(line, " vertex %f %f %f", &xx, &yy, &zz) == 3) {
                triVerts[vCount * 3 + 0] = xx;
                triVerts[vCount * 3 + 1] = yy;
                triVerts[vCount * 3 + 2] = zz;
                vCount++;
                if (vCount == 3) {
                    int va = addVertex(triVerts[0], triVerts[1], triVerts[2]);
                    int vb = addVertex(triVerts[3], triVerts[4], triVerts[5]);
                    int vc = addVertex(triVerts[6], triVerts[7], triVerts[8]);
                    t0.push_back(va);
                    t1.push_back(vb);
                    t2.push_back(vc);
                    vCount = 0;
                }
            }
        }
    }
    fclose(fp);

    int nVertices  = (int)vx.size();
    int nTriangles = (int)t0.size();
    printf("Deduplicated: %d vertices, %d triangles\n", nVertices, nTriangles);

    // Build adjacency, extract boundary edges
    std::vector<int> adj(nTriangles * 3, -1);
    std::vector<int> edgeVa, edgeVb;
    {
        std::unordered_map<uint64_t, int64_t> edgeToTri;
        edgeToTri.reserve(nTriangles * 3);
        int triVerts[3];
        for (int ii = 0; ii < nTriangles; ii++) {
            triVerts[0] = t0[ii]; triVerts[1] = t1[ii]; triVerts[2] = t2[ii];
            for (int ee = 0; ee < 3; ee++) {
                uint64_t key = edgeKey(triVerts[ee], triVerts[(ee + 1) % 3]);
                int64_t val = ((int64_t)ii << 2) | ee;
                auto [it, inserted] = edgeToTri.emplace(key, val);
                if (!inserted) {
                    int otherTri  = (int)(it->second >> 2);
                    int otherEdge = (int)(it->second & 3);
                    adj[ii * 3 + ee] = otherTri;
                    adj[otherTri * 3 + otherEdge] = ii;
                    it->second = -1;
                }
            }
        }
        for (auto &[key, val] : edgeToTri) {
            if (val >= 0) {
                int va = (int)(key & 0xFFFFFFFF);
                int vb = (int)(key >> 32);
                edgeVa.push_back(va);
                edgeVb.push_back(vb);
            }
        }
    }
    int nEdges = (int)edgeVa.size();
    printf("Boundary edges: %d\n", nEdges);

    // Flood fill connected components
    std::vector<int> comp(nTriangles, -1);
    std::vector<int> stack;
    int nComponents = 0;
    for (int ii = 0; ii < nTriangles; ii++) {
        if (comp[ii] >= 0) continue;
        int compId = nComponents++;
        stack.push_back(ii);
        comp[ii] = compId;
        while (!stack.empty()) {
            int tri = stack.back();
            stack.pop_back();
            for (int ee = 0; ee < 3; ee++) {
                int nb = adj[tri * 3 + ee];
                if (nb >= 0 && comp[nb] < 0) {
                    comp[nb] = compId;
                    stack.push_back(nb);
                }
            }
        }
    }
    printf("Connected components: %d\n", nComponents);

    // Write meshb (1-indexed)
    int64_t outMsh = GmfOpenMesh(outputFile, GmfWrite, 2, 3);
    if (!outMsh) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    GmfSetKwd(outMsh, GmfVertices, (int64_t)nVertices);
    for (int ii = 0; ii < nVertices; ii++)
        GmfSetLin(outMsh, GmfVertices, vx[ii], vy[ii], vz[ii], 0);

    GmfSetKwd(outMsh, GmfTriangles, (int64_t)nTriangles);
    for (int ii = 0; ii < nTriangles; ii++)
        GmfSetLin(outMsh, GmfTriangles, t0[ii] + 1, t1[ii] + 1, t2[ii] + 1, comp[ii] + 1);

    if (nEdges > 0) {
        GmfSetKwd(outMsh, GmfEdges, (int64_t)nEdges);
        for (int ii = 0; ii < nEdges; ii++)
            GmfSetLin(outMsh, GmfEdges, edgeVa[ii] + 1, edgeVb[ii] + 1, 1);
    }

    GmfCloseMesh(outMsh);
    printf("Wrote %s\n", outputFile);
    return 0;
}

// ---- meshb to STL ----

static int meshbToStl(const char *inputFile, const char *outputFile)
{
    int dim;
    int64_t fh = GmfOpenMesh(inputFile, GmfRead, &dim, &dim);
    if (!fh) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    int64_t nVerts  = GmfStatKwd(fh, GmfVertices);
    int64_t nTris   = GmfStatKwd(fh, GmfTriangles);
    int64_t nQuads  = GmfStatKwd(fh, GmfQuadrilaterals);
    int64_t nCorner = GmfStatKwd(fh, GmfCorners);
    int64_t nEdge   = GmfStatKwd(fh, GmfEdges);
    int64_t nTet    = GmfStatKwd(fh, GmfTetrahedra);
    int64_t nPyr    = GmfStatKwd(fh, GmfPyramids);
    int64_t nPrism  = GmfStatKwd(fh, GmfPrisms);
    int64_t nHex    = GmfStatKwd(fh, GmfHexahedra);

    if (nCorner > 0) printf("warning, element type not supported Corners\n");
    if (nEdge   > 0) printf("warning, element type not supported Edges\n");
    if (nTet    > 0) printf("warning, element type not supported Tetrahedra\n");
    if (nPyr    > 0) printf("warning, element type not supported Pyramids\n");
    if (nPrism  > 0) printf("warning, element type not supported Prisms\n");
    if (nHex    > 0) printf("warning, element type not supported Hexahedra\n");

    std::vector<double> vx(nVerts), vy(nVerts), vz(nVerts);
    GmfGotoKwd(fh, GmfVertices);
    for (int64_t ii = 0; ii < nVerts; ii++) {
        int ref;
        GmfGetLin(fh, GmfVertices, &vx[ii], &vy[ii], &vz[ii], &ref);
    }

    std::vector<int> t0(nTris), t1(nTris), t2(nTris), tRef(nTris);
    if (nTris > 0) {
        GmfGotoKwd(fh, GmfTriangles);
        for (int64_t ii = 0; ii < nTris; ii++)
            GmfGetLin(fh, GmfTriangles, &t0[ii], &t1[ii], &t2[ii], &tRef[ii]);
    }

    std::vector<int> q0(nQuads), q1(nQuads), q2(nQuads), q3(nQuads), qRef(nQuads);
    if (nQuads > 0) {
        GmfGotoKwd(fh, GmfQuadrilaterals);
        for (int64_t ii = 0; ii < nQuads; ii++)
            GmfGetLin(fh, GmfQuadrilaterals,
                      &q0[ii], &q1[ii], &q2[ii], &q3[ii], &qRef[ii]);
    }

    GmfCloseMesh(fh);

    printf("Read %ld vertices, %ld triangles, %ld quads\n",
           (long)nVerts, (long)nTris, (long)nQuads);

    // Write binary STL
    FILE *fp = fopen(outputFile, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    char header[80] = {};
    snprintf(header, 80, "Binary STL from %s", inputFile);
    fwrite(header, 1, 80, fp);
    uint32_t nTriOut = (uint32_t)(nTris + 2 * nQuads);
    fwrite(&nTriOut, 4, 1, fp);

    auto writeTri = [&](int ia, int ib, int ic, int ref) {
        double ax = vx[ib] - vx[ia], ay = vy[ib] - vy[ia], az = vz[ib] - vz[ia];
        double bx = vx[ic] - vx[ia], by = vy[ic] - vy[ia], bz = vz[ic] - vz[ia];
        float nx = (float)(ay * bz - az * by);
        float ny = (float)(az * bx - ax * bz);
        float nz = (float)(ax * by - ay * bx);
        fwrite(&nx, 4, 1, fp);
        fwrite(&ny, 4, 1, fp);
        fwrite(&nz, 4, 1, fp);
        float verts[9] = {
            (float)vx[ia], (float)vy[ia], (float)vz[ia],
            (float)vx[ib], (float)vy[ib], (float)vz[ib],
            (float)vx[ic], (float)vy[ic], (float)vz[ic]
        };
        fwrite(verts, sizeof(float), 9, fp);
        uint16_t attr = (uint16_t)ref;
        fwrite(&attr, 2, 1, fp);
    };

    for (int64_t ii = 0; ii < nTris; ii++)
        writeTri(t0[ii] - 1, t1[ii] - 1, t2[ii] - 1, tRef[ii]);

    for (int64_t ii = 0; ii < nQuads; ii++) {
        int ia = q0[ii] - 1, ib = q1[ii] - 1, ic = q2[ii] - 1, id = q3[ii] - 1;
        writeTri(ia, ib, ic, qRef[ii]);
        writeTri(ia, ic, id, qRef[ii]);
    }

    fclose(fp);
    printf("Wrote %s\n", outputFile);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 3) { printUsage(argv[0]); return argc < 2 ? 0 : 1; }

    const char *inputFile  = argv[1];
    const char *outputFile = argv[2];

    bool inputIsMeshb = endsWith(inputFile, ".meshb") || endsWith(inputFile, ".mesh");
    bool inputIsStl   = endsWith(inputFile, ".stl");

    if (inputIsMeshb)
        return meshbToStl(inputFile, outputFile);
    else if (inputIsStl)
        return stlToMeshb(inputFile, outputFile);
    else {
        fprintf(stderr, "Cannot determine direction from extension: %s\n", inputFile);
        fprintf(stderr, "Expected .stl or .meshb/.mesh\n");
        return 1;
    }
}
