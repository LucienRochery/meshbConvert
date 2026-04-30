#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <algorithm>

extern "C" {
#include <libmeshb7.h>
}

static uint64_t edgeKey(int va, int vb)
{
    uint32_t lo = (uint32_t)(va < vb ? va : vb);
    uint32_t hi = (uint32_t)(va < vb ? vb : va);
    return ((uint64_t)hi << 32) | lo;
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
    printf("Converts between UGRID and meshb surface mesh formats.\n");
    printf("Direction is detected from file extensions.\n");
    printf("\n");
    printf("  UGRID -> meshb : %s model.lb8.ugrid model.meshb\n", prog);
    printf("  meshb -> UGRID : %s model.meshb model.lb8.ugrid\n", prog);
    printf("\n");
    printf("UGRID to meshb:\n");
    printf("  - Reads lb8.ugrid (little-endian, 8-byte doubles, 4-byte ints)\n");
    printf("  - Extracts triangles, quads, boundary edges\n");
    printf("  - Assigns connected-component refs to triangles\n");
    printf("  - Volume elements are skipped\n");
    printf("\n");
    printf("meshb to UGRID:\n");
    printf("  - Writes lb8.ugrid from meshb triangles and quads\n");
    printf("  - Triangle/quad refs become surface IDs\n");
    printf("  - No volume elements are written\n");
}

// ---- UGRID to meshb ----

static int ugridToMeshb(const char *inputFile, const char *outputFile)
{
    FILE *fp = fopen(inputFile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    int32_t header[7];
    if (fread(header, sizeof(int32_t), 7, fp) != 7) {
        fprintf(stderr, "Failed to read UGRID header\n");
        fclose(fp);
        return 1;
    }
    int nNodes = header[0], nTri = header[1], nQuad = header[2];
    int nTet = header[3], nPyr = header[4], nPrism = header[5], nHex = header[6];

    printf("UGRID: %d nodes, %d tris, %d quads, %d tets, %d pyrs, %d prisms, %d hexes\n",
           nNodes, nTri, nQuad, nTet, nPyr, nPrism, nHex);

    std::vector<double> vx(nNodes), vy(nNodes), vz(nNodes);
    for (int ii = 0; ii < nNodes; ii++) {
        double xyz[3];
        fread(xyz, sizeof(double), 3, fp);
        vx[ii] = xyz[0]; vy[ii] = xyz[1]; vz[ii] = xyz[2];
    }

    std::vector<int32_t> triConn(nTri * 3);
    fread(triConn.data(), sizeof(int32_t), nTri * 3, fp);

    std::vector<int32_t> quadConn(nQuad * 4);
    if (nQuad > 0)
        fread(quadConn.data(), sizeof(int32_t), nQuad * 4, fp);

    // Skip volume elements
    fseek(fp, (long)nTet * 4 * sizeof(int32_t), SEEK_CUR);
    fseek(fp, (long)nPyr * 5 * sizeof(int32_t), SEEK_CUR);
    fseek(fp, (long)nPrism * 6 * sizeof(int32_t), SEEK_CUR);
    fseek(fp, (long)nHex * 8 * sizeof(int32_t), SEEK_CUR);

    std::vector<int32_t> surfIds(nTri + nQuad);
    fread(surfIds.data(), sizeof(int32_t), nTri + nQuad, fp);
    fclose(fp);

    // Build triangle adjacency, extract boundary edges
    std::vector<int> adj(nTri * 3, -1);
    std::vector<int> edgeVa, edgeVb;
    {
        std::unordered_map<uint64_t, int64_t> edgeToTri;
        edgeToTri.reserve(nTri * 3);
        for (int ii = 0; ii < nTri; ii++) {
            int tv[3] = {triConn[ii * 3], triConn[ii * 3 + 1], triConn[ii * 3 + 2]};
            for (int ee = 0; ee < 3; ee++) {
                uint64_t key = edgeKey(tv[ee], tv[(ee + 1) % 3]);
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
                edgeVa.push_back((int)(key & 0xFFFFFFFF));
                edgeVb.push_back((int)(key >> 32));
            }
        }
    }
    int nEdges = (int)edgeVa.size();
    printf("Boundary edges: %d\n", nEdges);

    // Flood fill connected components
    std::vector<int> comp(nTri, -1);
    std::vector<int> stack;
    int nComponents = 0;
    for (int ii = 0; ii < nTri; ii++) {
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

    // Write meshb (UGRID is already 1-indexed)
    int64_t outMsh = GmfOpenMesh(outputFile, GmfWrite, 2, 3);
    if (!outMsh) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    GmfSetKwd(outMsh, GmfVertices, (int64_t)nNodes);
    for (int ii = 0; ii < nNodes; ii++)
        GmfSetLin(outMsh, GmfVertices, vx[ii], vy[ii], vz[ii], 0);

    GmfSetKwd(outMsh, GmfTriangles, (int64_t)nTri);
    for (int ii = 0; ii < nTri; ii++)
        GmfSetLin(outMsh, GmfTriangles,
                  triConn[ii * 3], triConn[ii * 3 + 1], triConn[ii * 3 + 2],
                  comp[ii] + 1);

    if (nQuad > 0) {
        GmfSetKwd(outMsh, GmfQuadrilaterals, (int64_t)nQuad);
        for (int ii = 0; ii < nQuad; ii++)
            GmfSetLin(outMsh, GmfQuadrilaterals,
                      quadConn[ii * 4], quadConn[ii * 4 + 1],
                      quadConn[ii * 4 + 2], quadConn[ii * 4 + 3],
                      surfIds[nTri + ii]);
    }

    if (nEdges > 0) {
        GmfSetKwd(outMsh, GmfEdges, (int64_t)nEdges);
        for (int ii = 0; ii < nEdges; ii++)
            GmfSetLin(outMsh, GmfEdges, edgeVa[ii], edgeVb[ii], 0);
    }

    if (nTet > 0 || nPyr > 0 || nPrism > 0 || nHex > 0)
        printf("Note: volume elements (tet=%d pyr=%d prism=%d hex=%d) not written\n",
               nTet, nPyr, nPrism, nHex);

    GmfCloseMesh(outMsh);
    printf("Wrote %s\n", outputFile);
    return 0;
}

// ---- meshb to UGRID ----

static int meshbToUgrid(const char *inputFile, const char *outputFile)
{
    int dim;
    int64_t fh = GmfOpenMesh(inputFile, GmfRead, &dim, &dim);
    if (!fh) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    int64_t nNodes = GmfStatKwd(fh, GmfVertices);
    int64_t nTri   = GmfStatKwd(fh, GmfTriangles);
    int64_t nQuad  = GmfStatKwd(fh, GmfQuadrilaterals);

    std::vector<double> vx(nNodes), vy(nNodes), vz(nNodes);
    GmfGotoKwd(fh, GmfVertices);
    for (int64_t ii = 0; ii < nNodes; ii++) {
        int ref;
        GmfGetLin(fh, GmfVertices, &vx[ii], &vy[ii], &vz[ii], &ref);
    }

    std::vector<int32_t> triConn(nTri * 3), triRef(nTri);
    if (nTri > 0) {
        GmfGotoKwd(fh, GmfTriangles);
        for (int64_t ii = 0; ii < nTri; ii++) {
            int va, vb, vc, ref;
            GmfGetLin(fh, GmfTriangles, &va, &vb, &vc, &ref);
            triConn[ii * 3] = va; triConn[ii * 3 + 1] = vb; triConn[ii * 3 + 2] = vc;
            triRef[ii] = ref;
        }
    }

    std::vector<int32_t> quadConn(nQuad * 4), quadRef(nQuad);
    if (nQuad > 0) {
        GmfGotoKwd(fh, GmfQuadrilaterals);
        for (int64_t ii = 0; ii < nQuad; ii++) {
            int va, vb, vc, vd, ref;
            GmfGetLin(fh, GmfQuadrilaterals, &va, &vb, &vc, &vd, &ref);
            quadConn[ii * 4] = va; quadConn[ii * 4 + 1] = vb;
            quadConn[ii * 4 + 2] = vc; quadConn[ii * 4 + 3] = vd;
            quadRef[ii] = ref;
        }
    }

    GmfCloseMesh(fh);

    printf("meshb: %ld nodes, %ld tris, %ld quads\n", (long)nNodes, (long)nTri, (long)nQuad);

    // Compact refs to 1..N (sorted by original ref) so the ugrid surface IDs
    // and the .mapbc patch list are sequential and consistent.
    std::vector<int32_t> origRefs;
    origRefs.reserve(nTri + nQuad);
    for (int64_t ii = 0; ii < nTri; ii++) origRefs.push_back(triRef[ii]);
    for (int64_t ii = 0; ii < nQuad; ii++) origRefs.push_back(quadRef[ii]);
    std::sort(origRefs.begin(), origRefs.end());
    origRefs.erase(std::unique(origRefs.begin(), origRefs.end()), origRefs.end());

    std::unordered_map<int32_t, int32_t> refRemap;
    refRemap.reserve(origRefs.size() * 2);
    for (size_t ii = 0; ii < origRefs.size(); ii++)
        refRemap[origRefs[ii]] = (int32_t)(ii + 1);

    bool refsChanged = false;
    for (size_t ii = 0; ii < origRefs.size(); ii++)
        if (origRefs[ii] != (int32_t)(ii + 1)) { refsChanged = true; break; }
    if (refsChanged)
        printf("Compacting %d non-sequential refs to 1..%d\n",
               (int)origRefs.size(), (int)origRefs.size());

    for (int64_t ii = 0; ii < nTri; ii++) triRef[ii] = refRemap[triRef[ii]];
    for (int64_t ii = 0; ii < nQuad; ii++) quadRef[ii] = refRemap[quadRef[ii]];

    // Write lb8.ugrid: header, coords, tri conn, quad conn, (no vol elts), surface IDs
    FILE *fp = fopen(outputFile, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    int32_t header[7] = { (int32_t)nNodes, (int32_t)nTri, (int32_t)nQuad, 0, 0, 0, 0 };
    fwrite(header, sizeof(int32_t), 7, fp);

    for (int64_t ii = 0; ii < nNodes; ii++) {
        double xyz[3] = { vx[ii], vy[ii], vz[ii] };
        fwrite(xyz, sizeof(double), 3, fp);
    }

    fwrite(triConn.data(), sizeof(int32_t), nTri * 3, fp);
    if (nQuad > 0)
        fwrite(quadConn.data(), sizeof(int32_t), nQuad * 4, fp);

    // Surface IDs: tri refs then quad refs
    for (int64_t ii = 0; ii < nTri; ii++)
        fwrite(&triRef[ii], sizeof(int32_t), 1, fp);
    for (int64_t ii = 0; ii < nQuad; ii++)
        fwrite(&quadRef[ii], sizeof(int32_t), 1, fp);

    fclose(fp);
    printf("Wrote %s\n", outputFile);

    // Write companion .mapbc listing each unique surface ID
    char mapbcFile[4096];
    int lenOut = (int)strlen(outputFile);
    int stripLen = 0;
    if (endsWith(outputFile, ".lb8.ugrid")) stripLen = 10;
    else if (endsWith(outputFile, ".b8.ugrid")) stripLen = 9;
    else if (endsWith(outputFile, ".ugrid")) stripLen = 6;
    int baseLen = lenOut - stripLen;
    if (baseLen + 7 >= (int)sizeof(mapbcFile)) {
        fprintf(stderr, "Output path too long for .mapbc companion\n");
        return 1;
    }
    memcpy(mapbcFile, outputFile, baseLen);
    memcpy(mapbcFile + baseLen, ".mapbc", 7);

    FILE *fmap = fopen(mapbcFile, "w");
    if (!fmap) { fprintf(stderr, "Cannot open %s for writing\n", mapbcFile); return 1; }
    fprintf(fmap, "%d\n", (int)origRefs.size());
    for (size_t ii = 0; ii < origRefs.size(); ii++)
        fprintf(fmap, "%d 0 patch_%d\n", (int)(ii + 1), (int)(ii + 1));
    fclose(fmap);
    printf("Wrote %s (%d patches)\n", mapbcFile, (int)origRefs.size());

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 3) { printUsage(argv[0]); return argc < 2 ? 0 : 1; }

    const char *inputFile  = argv[1];
    const char *outputFile = argv[2];

    bool inputIsMeshb = endsWith(inputFile, ".meshb") || endsWith(inputFile, ".mesh");
    bool inputIsUgrid = endsWith(inputFile, ".ugrid");

    if (inputIsUgrid)
        return ugridToMeshb(inputFile, outputFile);
    else if (inputIsMeshb)
        return meshbToUgrid(inputFile, outputFile);
    else {
        fprintf(stderr, "Cannot determine direction from extension: %s\n", inputFile);
        fprintf(stderr, "Expected .ugrid or .meshb/.mesh\n");
        return 1;
    }
}
