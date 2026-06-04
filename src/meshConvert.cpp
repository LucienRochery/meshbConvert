#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cctype>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>

extern "C" {
#include <libmeshb7.h>
}

// ============================================================
// Common in-memory mesh (0-indexed; converted to/from 1-indexed at I/O)
// ============================================================

struct Mesh {
    std::vector<double> x, y, z;
    std::vector<int> triConn,   triRef;
    std::vector<int> quadConn,  quadRef;
    std::vector<int> tetConn,   tetRef;
    std::vector<int> pyrConn,   pyrRef;
    std::vector<int> prismConn, prismRef;
    std::vector<int> hexConn,   hexRef;
    std::vector<int> edgeConn,  edgeRef;

    int nVerts() const { return (int)x.size(); }
    int nTri()   const { return (int)(triConn.size()   / 3); }
    int nQuad()  const { return (int)(quadConn.size()  / 4); }
    int nTet()   const { return (int)(tetConn.size()   / 4); }
    int nPyr()   const { return (int)(pyrConn.size()   / 5); }
    int nPrism() const { return (int)(prismConn.size() / 6); }
    int nHex()   const { return (int)(hexConn.size()   / 8); }
    int nEdge()  const { return (int)(edgeConn.size()  / 2); }
};

static bool endsWith(const char *str, const char *suffix)
{
    size_t lenStr = strlen(str), lenSuf = strlen(suffix);
    if (lenSuf > lenStr) return false;
    return strcmp(str + lenStr - lenSuf, suffix) == 0;
}

static uint64_t edgeKey(int va, int vb)
{
    uint32_t lo = (uint32_t)(va < vb ? va : vb);
    uint32_t hi = (uint32_t)(va < vb ? vb : va);
    return ((uint64_t)hi << 32) | lo;
}

// ============================================================
// Derived: boundary edges and triangle connected components
// ============================================================

static void computeBoundaryEdges(const Mesh &mm, std::vector<int> &edgeConnOut)
{
    int nTri = mm.nTri();
    std::unordered_map<uint64_t, int> edgeCount;
    edgeCount.reserve(nTri * 3);
    for (int ii = 0; ii < nTri; ii++) {
        const int *tv = &mm.triConn[ii * 3];
        for (int ee = 0; ee < 3; ee++)
            edgeCount[edgeKey(tv[ee], tv[(ee + 1) % 3])]++;
    }
    for (auto &kv : edgeCount) {
        if (kv.second == 1) {
            edgeConnOut.push_back((int)(kv.first & 0xFFFFFFFF));
            edgeConnOut.push_back((int)(kv.first >> 32));
        }
    }
}

static int computeTriangleComponents(const Mesh &mm, std::vector<int> &compOut)
{
    int nTri = mm.nTri();
    std::vector<int> adj(nTri * 3, -1);
    std::unordered_map<uint64_t, int64_t> edgeToTri;
    edgeToTri.reserve(nTri * 3);
    for (int ii = 0; ii < nTri; ii++) {
        const int *tv = &mm.triConn[ii * 3];
        for (int ee = 0; ee < 3; ee++) {
            uint64_t kk = edgeKey(tv[ee], tv[(ee + 1) % 3]);
            int64_t val = ((int64_t)ii << 2) | ee;
            auto [it, inserted] = edgeToTri.emplace(kk, val);
            if (!inserted) {
                int otherTri  = (int)(it->second >> 2);
                int otherEdge = (int)(it->second & 3);
                adj[ii * 3 + ee] = otherTri;
                adj[otherTri * 3 + otherEdge] = ii;
                it->second = -1;
            }
        }
    }
    compOut.assign(nTri, -1);
    std::vector<int> stack;
    int nComp = 0;
    for (int ii = 0; ii < nTri; ii++) {
        if (compOut[ii] >= 0) continue;
        int cid = nComp++;
        compOut[ii] = cid;
        stack.push_back(ii);
        while (!stack.empty()) {
            int tt = stack.back(); stack.pop_back();
            for (int ee = 0; ee < 3; ee++) {
                int nb = adj[tt * 3 + ee];
                if (nb >= 0 && compOut[nb] < 0) {
                    compOut[nb] = cid;
                    stack.push_back(nb);
                }
            }
        }
    }
    return nComp;
}

// Fill in connected-component triangle refs (if all are zero) and boundary edges
// (if none). Runs unconditionally after read so every downstream writer sees
// the same enriched mesh.
static void fillMissingDerived(Mesh &mm)
{
    if (mm.nTri() == 0) return;
    if ((int)mm.triRef.size() != mm.nTri()) mm.triRef.assign(mm.nTri(), 0);

    bool allZero = true;
    for (int rr : mm.triRef) if (rr != 0) { allZero = false; break; }
    if (allZero) {
        std::vector<int> comp;
        int nComp = computeTriangleComponents(mm, comp);
        for (int ii = 0; ii < mm.nTri(); ii++) mm.triRef[ii] = comp[ii] + 1;
        printf("Assigned %d connected-component refs\n", nComp);
    }

    if (mm.nEdge() == 0) {
        computeBoundaryEdges(mm, mm.edgeConn);
        mm.edgeRef.assign(mm.nEdge(), 1);
        if (mm.nEdge() > 0) printf("Extracted %d boundary edges\n", mm.nEdge());
    }
}

// ============================================================
// STL
// ============================================================

struct Vec3Key {
    uint64_t bits[3];
    bool operator==(const Vec3Key &o) const
    { return bits[0]==o.bits[0] && bits[1]==o.bits[1] && bits[2]==o.bits[2]; }
};
struct Vec3KeyHash {
    size_t operator()(const Vec3Key &kk) const
    {
        size_t hh = 14695981039346656037ULL;
        for (int ii = 0; ii < 3; ii++) { hh ^= kk.bits[ii]; hh *= 1099511628211ULL; }
        return hh;
    }
};
static Vec3Key makeKey(float xx, float yy, float zz)
{
    Vec3Key kk{};
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
    long sz = ftell(fp);
    rewind(fp);
    return sz == (long)(80 + 4 + (long)nTri * 50);
}

static int readStl(const char *inputFile, Mesh &mm)
{
    FILE *fp = fopen(inputFile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    std::unordered_map<Vec3Key, int, Vec3KeyHash> vertexMap;
    auto addVertex = [&](float xx, float yy, float zz) -> int {
        Vec3Key kk = makeKey(xx, yy, zz);
        auto [it, inserted] = vertexMap.emplace(kk, (int)mm.x.size());
        if (inserted) {
            mm.x.push_back((double)xx); mm.y.push_back((double)yy); mm.z.push_back((double)zz);
        }
        return it->second;
    };

    if (isBinaryStl(fp)) {
        fseek(fp, 80, SEEK_SET);
        uint32_t nTri;
        fread(&nTri, 4, 1, fp);
        printf("Reading binary STL: %u triangles\n", nTri);
        for (uint32_t ii = 0; ii < nTri; ii++) {
            float buf[12]; uint16_t attr;
            fread(buf, sizeof(float), 12, fp);
            fread(&attr, 2, 1, fp);
            int va = addVertex(buf[3], buf[4], buf[5]);
            int vb = addVertex(buf[6], buf[7], buf[8]);
            int vc = addVertex(buf[9], buf[10], buf[11]);
            mm.triConn.push_back(va); mm.triConn.push_back(vb); mm.triConn.push_back(vc);
        }
    } else {
        fseek(fp, 0, SEEK_SET);
        char line[256]; float tv[9]; int vCount = 0;
        while (fgets(line, sizeof(line), fp)) {
            float xx, yy, zz;
            if (sscanf(line, " vertex %f %f %f", &xx, &yy, &zz) == 3) {
                tv[vCount*3+0]=xx; tv[vCount*3+1]=yy; tv[vCount*3+2]=zz;
                vCount++;
                if (vCount == 3) {
                    int va = addVertex(tv[0], tv[1], tv[2]);
                    int vb = addVertex(tv[3], tv[4], tv[5]);
                    int vc = addVertex(tv[6], tv[7], tv[8]);
                    mm.triConn.push_back(va); mm.triConn.push_back(vb); mm.triConn.push_back(vc);
                    vCount = 0;
                }
            }
        }
    }
    fclose(fp);

    mm.triRef.assign(mm.nTri(), 0);
    printf("STL: %d verts, %d triangles (deduplicated)\n", mm.nVerts(), mm.nTri());
    return 0;
}

static int writeStl(const char *outputFile, const Mesh &mm)
{
    FILE *fp = fopen(outputFile, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    char header[80] = {};
    snprintf(header, 80, "Binary STL");
    fwrite(header, 1, 80, fp);

    uint32_t nOut = (uint32_t)(mm.nTri() + 2 * mm.nQuad());
    fwrite(&nOut, 4, 1, fp);

    auto emitTri = [&](int ia, int ib, int ic, uint16_t attr) {
        double ax = mm.x[ib] - mm.x[ia], ay = mm.y[ib] - mm.y[ia], az = mm.z[ib] - mm.z[ia];
        double bx = mm.x[ic] - mm.x[ia], by = mm.y[ic] - mm.y[ia], bz = mm.z[ic] - mm.z[ia];
        float nn[3] = {
            (float)(ay * bz - az * by),
            (float)(az * bx - ax * bz),
            (float)(ax * by - ay * bx)
        };
        fwrite(nn, 4, 3, fp);
        float verts[9] = {
            (float)mm.x[ia], (float)mm.y[ia], (float)mm.z[ia],
            (float)mm.x[ib], (float)mm.y[ib], (float)mm.z[ib],
            (float)mm.x[ic], (float)mm.y[ic], (float)mm.z[ic]
        };
        fwrite(verts, 4, 9, fp);
        fwrite(&attr, 2, 1, fp);
    };

    for (int ii = 0; ii < mm.nTri(); ii++) {
        uint16_t attr = (uint16_t)(mm.triRef.empty() ? 0 : mm.triRef[ii]);
        emitTri(mm.triConn[ii*3], mm.triConn[ii*3+1], mm.triConn[ii*3+2], attr);
    }
    for (int ii = 0; ii < mm.nQuad(); ii++) {
        uint16_t attr = (uint16_t)(mm.quadRef.empty() ? 0 : mm.quadRef[ii]);
        int qv[4] = { mm.quadConn[ii*4], mm.quadConn[ii*4+1], mm.quadConn[ii*4+2], mm.quadConn[ii*4+3] };
        emitTri(qv[0], qv[1], qv[2], attr);
        emitTri(qv[0], qv[2], qv[3], attr);
    }
    if (mm.nTet() + mm.nPyr() + mm.nPrism() + mm.nHex() > 0)
        printf("Note: %d volume cells not written to STL\n",
               mm.nTet() + mm.nPyr() + mm.nPrism() + mm.nHex());

    fclose(fp);
    printf("Wrote %s (%u triangles)\n", outputFile, nOut);
    return 0;
}

// ============================================================
// meshb
// ============================================================

static int readMeshb(const char *inputFile, Mesh &mm)
{
    int ver, dim;
    int64_t fh = GmfOpenMesh(inputFile, GmfRead, &ver, &dim);
    if (!fh) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    int64_t nV  = GmfStatKwd(fh, GmfVertices);
    int64_t nT  = GmfStatKwd(fh, GmfTriangles);
    int64_t nQ  = GmfStatKwd(fh, GmfQuadrilaterals);
    int64_t nTe = GmfStatKwd(fh, GmfTetrahedra);
    int64_t nPy = GmfStatKwd(fh, GmfPyramids);
    int64_t nPr = GmfStatKwd(fh, GmfPrisms);
    int64_t nHe = GmfStatKwd(fh, GmfHexahedra);
    int64_t nE  = GmfStatKwd(fh, GmfEdges);

    mm.x.resize(nV); mm.y.resize(nV); mm.z.resize(nV);
    GmfGotoKwd(fh, GmfVertices);
    for (int64_t ii = 0; ii < nV; ii++) {
        int ref;
        GmfGetLin(fh, GmfVertices, &mm.x[ii], &mm.y[ii], &mm.z[ii], &ref);
    }

    if (nT > 0) {
        mm.triConn.resize(nT * 3); mm.triRef.resize(nT);
        GmfGotoKwd(fh, GmfTriangles);
        for (int64_t ii = 0; ii < nT; ii++) {
            int v[3], rr;
            GmfGetLin(fh, GmfTriangles, &v[0], &v[1], &v[2], &rr);
            mm.triConn[ii*3]=v[0]-1; mm.triConn[ii*3+1]=v[1]-1; mm.triConn[ii*3+2]=v[2]-1;
            mm.triRef[ii] = rr;
        }
    }
    if (nQ > 0) {
        mm.quadConn.resize(nQ * 4); mm.quadRef.resize(nQ);
        GmfGotoKwd(fh, GmfQuadrilaterals);
        for (int64_t ii = 0; ii < nQ; ii++) {
            int v[4], rr;
            GmfGetLin(fh, GmfQuadrilaterals, &v[0], &v[1], &v[2], &v[3], &rr);
            for (int kk = 0; kk < 4; kk++) mm.quadConn[ii*4+kk] = v[kk] - 1;
            mm.quadRef[ii] = rr;
        }
    }
    if (nTe > 0) {
        mm.tetConn.resize(nTe * 4); mm.tetRef.resize(nTe);
        GmfGotoKwd(fh, GmfTetrahedra);
        for (int64_t ii = 0; ii < nTe; ii++) {
            int v[4], rr;
            GmfGetLin(fh, GmfTetrahedra, &v[0], &v[1], &v[2], &v[3], &rr);
            for (int kk = 0; kk < 4; kk++) mm.tetConn[ii*4+kk] = v[kk] - 1;
            mm.tetRef[ii] = rr;
        }
    }
    if (nPy > 0) {
        mm.pyrConn.resize(nPy * 5); mm.pyrRef.resize(nPy);
        GmfGotoKwd(fh, GmfPyramids);
        for (int64_t ii = 0; ii < nPy; ii++) {
            int v[5], rr;
            GmfGetLin(fh, GmfPyramids, &v[0], &v[1], &v[2], &v[3], &v[4], &rr);
            for (int kk = 0; kk < 5; kk++) mm.pyrConn[ii*5+kk] = v[kk] - 1;
            mm.pyrRef[ii] = rr;
        }
    }
    if (nPr > 0) {
        mm.prismConn.resize(nPr * 6); mm.prismRef.resize(nPr);
        GmfGotoKwd(fh, GmfPrisms);
        for (int64_t ii = 0; ii < nPr; ii++) {
            int v[6], rr;
            GmfGetLin(fh, GmfPrisms, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &rr);
            for (int kk = 0; kk < 6; kk++) mm.prismConn[ii*6+kk] = v[kk] - 1;
            mm.prismRef[ii] = rr;
        }
    }
    if (nHe > 0) {
        mm.hexConn.resize(nHe * 8); mm.hexRef.resize(nHe);
        GmfGotoKwd(fh, GmfHexahedra);
        for (int64_t ii = 0; ii < nHe; ii++) {
            int v[8], rr;
            GmfGetLin(fh, GmfHexahedra,
                      &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &rr);
            for (int kk = 0; kk < 8; kk++) mm.hexConn[ii*8+kk] = v[kk] - 1;
            mm.hexRef[ii] = rr;
        }
    }
    if (nE > 0) {
        mm.edgeConn.resize(nE * 2); mm.edgeRef.resize(nE);
        GmfGotoKwd(fh, GmfEdges);
        for (int64_t ii = 0; ii < nE; ii++) {
            int v[2], rr;
            GmfGetLin(fh, GmfEdges, &v[0], &v[1], &rr);
            mm.edgeConn[ii*2]=v[0]-1; mm.edgeConn[ii*2+1]=v[1]-1;
            mm.edgeRef[ii] = rr;
        }
    }
    GmfCloseMesh(fh);

    printf("meshb: %d verts, %d tri, %d quad, %d tet, %d pyr, %d prism, %d hex, %d edges\n",
           mm.nVerts(), mm.nTri(), mm.nQuad(), mm.nTet(),
           mm.nPyr(), mm.nPrism(), mm.nHex(), mm.nEdge());
    return 0;
}

static int writeMeshb(const char *outputFile, const Mesh &mm)
{
    int64_t fh = GmfOpenMesh(outputFile, GmfWrite, 2, 3);
    if (!fh) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    GmfSetKwd(fh, GmfVertices, (int64_t)mm.nVerts());
    for (int ii = 0; ii < mm.nVerts(); ii++)
        GmfSetLin(fh, GmfVertices, mm.x[ii], mm.y[ii], mm.z[ii], 0);

    if (mm.nTri() > 0) {
        GmfSetKwd(fh, GmfTriangles, (int64_t)mm.nTri());
        for (int ii = 0; ii < mm.nTri(); ii++)
            GmfSetLin(fh, GmfTriangles,
                      mm.triConn[ii*3]+1, mm.triConn[ii*3+1]+1, mm.triConn[ii*3+2]+1,
                      mm.triRef[ii]);
    }
    if (mm.nQuad() > 0) {
        GmfSetKwd(fh, GmfQuadrilaterals, (int64_t)mm.nQuad());
        for (int ii = 0; ii < mm.nQuad(); ii++)
            GmfSetLin(fh, GmfQuadrilaterals,
                      mm.quadConn[ii*4]+1, mm.quadConn[ii*4+1]+1,
                      mm.quadConn[ii*4+2]+1, mm.quadConn[ii*4+3]+1,
                      mm.quadRef.empty() ? 0 : mm.quadRef[ii]);
    }
    if (mm.nTet() > 0) {
        GmfSetKwd(fh, GmfTetrahedra, (int64_t)mm.nTet());
        for (int ii = 0; ii < mm.nTet(); ii++)
            GmfSetLin(fh, GmfTetrahedra,
                      mm.tetConn[ii*4]+1, mm.tetConn[ii*4+1]+1,
                      mm.tetConn[ii*4+2]+1, mm.tetConn[ii*4+3]+1,
                      mm.tetRef.empty() ? 0 : mm.tetRef[ii]);
    }
    if (mm.nPyr() > 0) {
        GmfSetKwd(fh, GmfPyramids, (int64_t)mm.nPyr());
        for (int ii = 0; ii < mm.nPyr(); ii++)
            GmfSetLin(fh, GmfPyramids,
                      mm.pyrConn[ii*5]+1, mm.pyrConn[ii*5+1]+1,
                      mm.pyrConn[ii*5+2]+1, mm.pyrConn[ii*5+3]+1,
                      mm.pyrConn[ii*5+4]+1,
                      mm.pyrRef.empty() ? 0 : mm.pyrRef[ii]);
    }
    if (mm.nPrism() > 0) {
        GmfSetKwd(fh, GmfPrisms, (int64_t)mm.nPrism());
        for (int ii = 0; ii < mm.nPrism(); ii++)
            GmfSetLin(fh, GmfPrisms,
                      mm.prismConn[ii*6]+1, mm.prismConn[ii*6+1]+1,
                      mm.prismConn[ii*6+2]+1, mm.prismConn[ii*6+3]+1,
                      mm.prismConn[ii*6+4]+1, mm.prismConn[ii*6+5]+1,
                      mm.prismRef.empty() ? 0 : mm.prismRef[ii]);
    }
    if (mm.nHex() > 0) {
        GmfSetKwd(fh, GmfHexahedra, (int64_t)mm.nHex());
        for (int ii = 0; ii < mm.nHex(); ii++)
            GmfSetLin(fh, GmfHexahedra,
                      mm.hexConn[ii*8]+1, mm.hexConn[ii*8+1]+1,
                      mm.hexConn[ii*8+2]+1, mm.hexConn[ii*8+3]+1,
                      mm.hexConn[ii*8+4]+1, mm.hexConn[ii*8+5]+1,
                      mm.hexConn[ii*8+6]+1, mm.hexConn[ii*8+7]+1,
                      mm.hexRef.empty() ? 0 : mm.hexRef[ii]);
    }
    if (mm.nEdge() > 0) {
        GmfSetKwd(fh, GmfEdges, (int64_t)mm.nEdge());
        for (int ii = 0; ii < mm.nEdge(); ii++)
            GmfSetLin(fh, GmfEdges,
                      mm.edgeConn[ii*2]+1, mm.edgeConn[ii*2+1]+1,
                      mm.edgeRef.empty() ? 1 : mm.edgeRef[ii]);
    }

    GmfCloseMesh(fh);
    printf("Wrote %s\n", outputFile);
    return 0;
}

// ============================================================
// UGRID (lb8.ugrid: little-endian, 8-byte doubles, 4-byte ints)
// ============================================================

static int readUgrid(const char *inputFile, Mesh &mm)
{
    FILE *fp = fopen(inputFile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }

    int32_t header[7];
    if (fread(header, sizeof(int32_t), 7, fp) != 7) {
        fprintf(stderr, "Failed to read UGRID header\n");
        fclose(fp); return 1;
    }
    int nNodes = header[0], nTri = header[1], nQuad = header[2];
    int nTet = header[3], nPyr = header[4], nPrism = header[5], nHex = header[6];
    printf("UGRID: %d nodes, %d tris, %d quads, %d tets, %d pyrs, %d prisms, %d hexes\n",
           nNodes, nTri, nQuad, nTet, nPyr, nPrism, nHex);

    mm.x.resize(nNodes); mm.y.resize(nNodes); mm.z.resize(nNodes);
    for (int ii = 0; ii < nNodes; ii++) {
        double xyz[3]; fread(xyz, sizeof(double), 3, fp);
        mm.x[ii] = xyz[0]; mm.y[ii] = xyz[1]; mm.z[ii] = xyz[2];
    }

    auto readConn = [&](int nElt, int stride, std::vector<int> &outConn) {
        if (nElt <= 0) return;
        std::vector<int32_t> tmp(nElt * stride);
        fread(tmp.data(), sizeof(int32_t), nElt * stride, fp);
        outConn.resize(nElt * stride);
        for (size_t ii = 0; ii < tmp.size(); ii++) outConn[ii] = tmp[ii] - 1;
    };
    readConn(nTri,   3, mm.triConn);
    readConn(nQuad,  4, mm.quadConn);
    readConn(nTet,   4, mm.tetConn);
    readConn(nPyr,   5, mm.pyrConn);
    readConn(nPrism, 6, mm.prismConn);
    readConn(nHex,   8, mm.hexConn);

    std::vector<int32_t> surfIds(nTri + nQuad);
    fread(surfIds.data(), sizeof(int32_t), nTri + nQuad, fp);
    fclose(fp);

    mm.triRef.assign(nTri, 0);
    mm.quadRef.assign(nQuad, 0);
    for (int ii = 0; ii < nTri;  ii++) mm.triRef[ii]  = surfIds[ii];
    for (int ii = 0; ii < nQuad; ii++) mm.quadRef[ii] = surfIds[nTri + ii];
    mm.tetRef.assign(nTet, 0);
    mm.pyrRef.assign(nPyr, 0);
    mm.prismRef.assign(nPrism, 0);
    mm.hexRef.assign(nHex, 0);
    return 0;
}

static int writeUgrid(const char *outputFile, const Mesh &mm)
{
    // Compact surface refs (tri+quad) to 1..N for the .mapbc companion
    std::vector<int32_t> origRefs;
    origRefs.reserve(mm.nTri() + mm.nQuad());
    for (int ii = 0; ii < mm.nTri();  ii++) origRefs.push_back(mm.triRef.empty()  ? 0 : mm.triRef[ii]);
    for (int ii = 0; ii < mm.nQuad(); ii++) origRefs.push_back(mm.quadRef.empty() ? 0 : mm.quadRef[ii]);
    std::sort(origRefs.begin(), origRefs.end());
    origRefs.erase(std::unique(origRefs.begin(), origRefs.end()), origRefs.end());

    std::unordered_map<int32_t, int32_t> remap;
    remap.reserve(origRefs.size() * 2);
    for (size_t ii = 0; ii < origRefs.size(); ii++) remap[origRefs[ii]] = (int32_t)(ii + 1);

    FILE *fp = fopen(outputFile, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    int32_t hdr[7] = {
        (int32_t)mm.nVerts(), (int32_t)mm.nTri(), (int32_t)mm.nQuad(),
        (int32_t)mm.nTet(),   (int32_t)mm.nPyr(), (int32_t)mm.nPrism(), (int32_t)mm.nHex()
    };
    fwrite(hdr, sizeof(int32_t), 7, fp);

    for (int ii = 0; ii < mm.nVerts(); ii++) {
        double xyz[3] = { mm.x[ii], mm.y[ii], mm.z[ii] };
        fwrite(xyz, sizeof(double), 3, fp);
    }

    auto dumpConn = [&](const std::vector<int> &conn) {
        if (conn.empty()) return;
        std::vector<int32_t> oneBased(conn.size());
        for (size_t ii = 0; ii < conn.size(); ii++) oneBased[ii] = conn[ii] + 1;
        fwrite(oneBased.data(), sizeof(int32_t), oneBased.size(), fp);
    };
    dumpConn(mm.triConn);
    dumpConn(mm.quadConn);
    dumpConn(mm.tetConn);
    dumpConn(mm.pyrConn);
    dumpConn(mm.prismConn);
    dumpConn(mm.hexConn);

    for (int ii = 0; ii < mm.nTri(); ii++) {
        int32_t rr = remap[mm.triRef.empty() ? 0 : mm.triRef[ii]];
        fwrite(&rr, sizeof(int32_t), 1, fp);
    }
    for (int ii = 0; ii < mm.nQuad(); ii++) {
        int32_t rr = remap[mm.quadRef.empty() ? 0 : mm.quadRef[ii]];
        fwrite(&rr, sizeof(int32_t), 1, fp);
    }
    fclose(fp);
    printf("Wrote %s\n", outputFile);

    // Companion .mapbc
    char mapbcFile[4096];
    size_t lenOut = strlen(outputFile);
    size_t stripLen = 0;
    if      (endsWith(outputFile, ".lb8.ugrid")) stripLen = 10;
    else if (endsWith(outputFile, ".b8.ugrid"))  stripLen = 9;
    else if (endsWith(outputFile, ".ugrid"))     stripLen = 6;
    size_t baseLen = lenOut - stripLen;
    if (baseLen + 7 >= sizeof(mapbcFile)) {
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

// ============================================================
// VTU (VTK UnstructuredGrid XML)
// ============================================================
// Supports format="ascii" and format="binary" (inline base64, uncompressed).
// Rejects format="appended" and any compressor.

struct XmlAttrs {
    std::vector<std::pair<std::string, std::string>> kv;
    std::string get(const char *key, const char *def = "") const
    {
        for (auto &pp : kv) if (pp.first == key) return pp.second;
        return std::string(def);
    }
};

// Scans the next opening / closing / self-closing tag at or after pos.
// Updates pos to just after '>'. Returns false at EOF.
static bool nextXmlTag(const std::string &buf, size_t &pos,
                       std::string &name, XmlAttrs &attrs,
                       bool &isClose, bool &isSelfClosing)
{
    while (pos < buf.size()) {
        size_t lt = buf.find('<', pos);
        if (lt == std::string::npos) return false;

        // Comment <!-- ... -->
        if (lt + 3 < buf.size() && buf[lt+1]=='!' && buf[lt+2]=='-' && buf[lt+3]=='-') {
            size_t end = buf.find("-->", lt);
            if (end == std::string::npos) return false;
            pos = end + 3;
            continue;
        }
        // CDATA <![CDATA[ ... ]]>
        if (lt + 8 < buf.size() && buf.compare(lt+1, 8, "![CDATA[") == 0) {
            size_t end = buf.find("]]>", lt);
            if (end == std::string::npos) return false;
            pos = end + 3;
            continue;
        }
        // Processing instruction <? ... ?>
        if (buf[lt+1] == '?') {
            size_t end = buf.find("?>", lt);
            if (end == std::string::npos) return false;
            pos = end + 2;
            continue;
        }

        size_t gt = buf.find('>', lt);
        if (gt == std::string::npos) return false;

        size_t pp = lt + 1;
        isClose = false;
        if (buf[pp] == '/') { isClose = true; pp++; }

        size_t nameStart = pp;
        while (pp < gt && (isalnum((unsigned char)buf[pp]) || buf[pp]==':' || buf[pp]=='_' || buf[pp]=='-'))
            pp++;
        name.assign(buf.data() + nameStart, pp - nameStart);

        attrs.kv.clear();
        while (pp < gt) {
            while (pp < gt && isspace((unsigned char)buf[pp])) pp++;
            if (pp >= gt || buf[pp] == '/') break;
            size_t kStart = pp;
            while (pp < gt && buf[pp] != '=' && !isspace((unsigned char)buf[pp])) pp++;
            std::string kk(buf.data() + kStart, pp - kStart);
            while (pp < gt && (isspace((unsigned char)buf[pp]) || buf[pp]=='=')) pp++;
            if (pp >= gt) break;
            char quote = buf[pp];
            if (quote != '"' && quote != '\'') break;
            pp++;
            size_t vStart = pp;
            while (pp < gt && buf[pp] != quote) pp++;
            std::string vv(buf.data() + vStart, pp - vStart);
            attrs.kv.push_back({kk, vv});
            if (pp < gt) pp++;
        }
        isSelfClosing = (gt > lt + 1 && buf[gt-1] == '/');
        pos = gt + 1;
        return true;
    }
    return false;
}

static void base64Decode(const std::string &ss, std::vector<uint8_t> &out)
{
    static int8_t dec[256];
    static bool init = false;
    if (!init) {
        for (int ii = 0; ii < 256; ii++) dec[ii] = -1;
        const char *abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        for (int ii = 0; ii < 64; ii++) dec[(unsigned char)abc[ii]] = (int8_t)ii;
        init = true;
    }
    out.clear();
    out.reserve(ss.size() * 3 / 4);
    // VTK concatenates two base64 chunks per DataArray (header then data),
    // each with its own '=' padding. Reset the bit accumulator on '=' so the
    // next chunk starts byte-aligned.
    uint32_t acc = 0; int bits = 0;
    for (char cc : ss) {
        if (cc == '=') { acc = 0; bits = 0; continue; }
        if (isspace((unsigned char)cc)) continue;
        int8_t vv = dec[(unsigned char)cc];
        if (vv < 0) continue;
        acc = (acc << 6) | (uint32_t)vv;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back((uint8_t)((acc >> bits) & 0xff));
        }
    }
}

template <typename T>
static void appendTyped(const uint8_t *data, size_t nBytes, const std::string &type,
                        std::vector<T> &out)
{
    if (type == "Float32") {
        size_t nn = nBytes / 4; const float *pp = (const float *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "Float64") {
        size_t nn = nBytes / 8; const double *pp = (const double *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "Int8") {
        for (size_t ii = 0; ii < nBytes; ii++) out.push_back((T)(int8_t)data[ii]);
    } else if (type == "UInt8") {
        for (size_t ii = 0; ii < nBytes; ii++) out.push_back((T)data[ii]);
    } else if (type == "Int16") {
        size_t nn = nBytes / 2; const int16_t *pp = (const int16_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "UInt16") {
        size_t nn = nBytes / 2; const uint16_t *pp = (const uint16_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "Int32") {
        size_t nn = nBytes / 4; const int32_t *pp = (const int32_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "UInt32") {
        size_t nn = nBytes / 4; const uint32_t *pp = (const uint32_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "Int64") {
        size_t nn = nBytes / 8; const int64_t *pp = (const int64_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else if (type == "UInt64") {
        size_t nn = nBytes / 8; const uint64_t *pp = (const uint64_t *)data;
        for (size_t ii = 0; ii < nn; ii++) out.push_back((T)pp[ii]);
    } else {
        fprintf(stderr, "VTU: unknown DataArray type '%s'\n", type.c_str());
    }
}

template <typename T>
static int decodeArray(const std::string &content, const std::string &format,
                       const std::string &type, int headerBytes,
                       std::vector<T> &out)
{
    out.clear();
    if (format == "ascii") {
        const char *pp = content.c_str();
        while (*pp) {
            while (*pp && isspace((unsigned char)*pp)) pp++;
            if (!*pp) break;
            char *endp;
            double val = strtod(pp, &endp);
            if (endp == pp) break;
            out.push_back((T)val);
            pp = endp;
        }
    } else if (format == "binary") {
        std::vector<uint8_t> raw;
        base64Decode(content, raw);
        if ((int)raw.size() < headerBytes) {
            fprintf(stderr, "VTU: binary block too small\n"); return 1;
        }
        uint64_t payloadBytes;
        if (headerBytes == 4) {
            uint32_t hh; memcpy(&hh, raw.data(), 4); payloadBytes = hh;
        } else {
            uint64_t hh; memcpy(&hh, raw.data(), 8); payloadBytes = hh;
        }
        if (raw.size() - headerBytes < payloadBytes) {
            fprintf(stderr, "VTU: binary payload short (%zu vs %llu)\n",
                    raw.size() - headerBytes, (unsigned long long)payloadBytes);
            return 1;
        }
        appendTyped(raw.data() + headerBytes, (size_t)payloadBytes, type, out);
    } else {
        fprintf(stderr, "VTU: format '%s' not supported\n", format.c_str());
        return 1;
    }
    return 0;
}

static int readVtu(const char *inputFile, Mesh &mm)
{
    FILE *fp = fopen(inputFile, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", inputFile); return 1; }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    rewind(fp);
    std::string buf(sz, '\0');
    fread(buf.data(), 1, sz, fp);
    fclose(fp);

    int headerBytes = 4;
    std::vector<double> pts;
    std::vector<int64_t> conn, offsets;
    std::vector<int> cellTypes;
    std::vector<int> cellRefs;
    std::string cellRefsName;
    int nPiecesSeen = 0;

    std::vector<std::string> stack;
    size_t pos = 0;
    std::string name; XmlAttrs attrs; bool isClose, isSelfClosing;

    while (nextXmlTag(buf, pos, name, attrs, isClose, isSelfClosing)) {
        if (isClose) {
            if (!stack.empty() && stack.back() == name) stack.pop_back();
            continue;
        }

        if (name == "VTKFile") {
            std::string ht = attrs.get("header_type", "UInt32");
            headerBytes = (ht == "UInt64") ? 8 : 4;
            std::string compr = attrs.get("compressor", "");
            if (!compr.empty()) {
                fprintf(stderr, "VTU: compressor '%s' not supported\n", compr.c_str());
                return 1;
            }
            std::string bo = attrs.get("byte_order", "LittleEndian");
            if (bo != "LittleEndian")
                fprintf(stderr, "VTU: byte_order='%s' not supported, assuming host endian\n", bo.c_str());
        }
        if (name == "Piece") nPiecesSeen++;

        if (name == "DataArray") {
            std::string daType   = attrs.get("type");
            std::string daName   = attrs.get("Name");
            std::string daFormat = attrs.get("format", "ascii");
            int daNComp = atoi(attrs.get("NumberOfComponents", "1").c_str());

            if (daFormat == "appended") {
                fprintf(stderr, "VTU: appended format not supported (DataArray '%s')\n", daName.c_str());
                return 1;
            }

            std::string content;
            if (!isSelfClosing) {
                size_t closeStart = buf.find("</DataArray", pos);
                if (closeStart == std::string::npos) {
                    fprintf(stderr, "VTU: unterminated DataArray '%s'\n", daName.c_str());
                    return 1;
                }
                content.assign(buf.data() + pos, closeStart - pos);
                size_t closeEnd = buf.find('>', closeStart);
                pos = closeEnd + 1;
            }
            if (nPiecesSeen > 1) continue;  // only first piece for now

            std::string parent = stack.empty() ? std::string() : stack.back();
            if (parent == "Points") {
                if (daNComp != 3) {
                    fprintf(stderr, "VTU: Points NumberOfComponents=%d, expected 3\n", daNComp);
                    return 1;
                }
                if (decodeArray(content, daFormat, daType, headerBytes, pts)) return 1;
            } else if (parent == "Cells") {
                if (daName == "connectivity") {
                    if (decodeArray(content, daFormat, daType, headerBytes, conn)) return 1;
                } else if (daName == "offsets") {
                    if (decodeArray(content, daFormat, daType, headerBytes, offsets)) return 1;
                } else if (daName == "types") {
                    std::vector<int64_t> tmp;
                    if (decodeArray(content, daFormat, daType, headerBytes, tmp)) return 1;
                    cellTypes.assign(tmp.begin(), tmp.end());
                }
            } else if (parent == "CellData" && cellRefs.empty()) {
                bool isInt = (daType.find("Int") != std::string::npos);
                if (isInt && daNComp == 1) {
                    std::vector<int64_t> tmp;
                    if (decodeArray(content, daFormat, daType, headerBytes, tmp)) return 1;
                    cellRefs.assign(tmp.begin(), tmp.end());
                    cellRefsName = daName;
                }
            }
        } else if (!isSelfClosing) {
            stack.push_back(name);
        }
    }

    if (nPiecesSeen > 1)
        fprintf(stderr, "VTU: %d Pieces in file; only first piece is read\n", nPiecesSeen);

    int nVerts = (int)(pts.size() / 3);
    mm.x.resize(nVerts); mm.y.resize(nVerts); mm.z.resize(nVerts);
    for (int ii = 0; ii < nVerts; ii++) {
        mm.x[ii] = pts[ii*3]; mm.y[ii] = pts[ii*3+1]; mm.z[ii] = pts[ii*3+2];
    }

    int nCells = (int)cellTypes.size();
    if ((int)offsets.size() != nCells) {
        fprintf(stderr, "VTU: offsets/types size mismatch (%d vs %d)\n",
                (int)offsets.size(), nCells);
        return 1;
    }
    bool hasRefs = (int)cellRefs.size() == nCells;
    if (!cellRefs.empty() && !hasRefs) {
        fprintf(stderr, "VTU: CellData '%s' size %d != nCells %d; ignoring\n",
                cellRefsName.c_str(), (int)cellRefs.size(), nCells);
    }

    int64_t prev = 0;
    int nUnhandled = 0;
    for (int ii = 0; ii < nCells; ii++) {
        int64_t cur = offsets[ii];
        int nn = (int)(cur - prev);
        const int64_t *vv = &conn[prev];
        int rr = hasRefs ? cellRefs[ii] : 0;
        switch (cellTypes[ii]) {
            case 5: // VTK_TRIANGLE
                if (nn == 3) {
                    mm.triConn.push_back((int)vv[0]);
                    mm.triConn.push_back((int)vv[1]);
                    mm.triConn.push_back((int)vv[2]);
                    mm.triRef.push_back(rr);
                } else nUnhandled++;
                break;
            case 9: // VTK_QUAD
                if (nn == 4) {
                    for (int kk = 0; kk < 4; kk++) mm.quadConn.push_back((int)vv[kk]);
                    mm.quadRef.push_back(rr);
                } else nUnhandled++;
                break;
            case 10: // VTK_TETRA
                if (nn == 4) {
                    for (int kk = 0; kk < 4; kk++) mm.tetConn.push_back((int)vv[kk]);
                    mm.tetRef.push_back(rr);
                } else nUnhandled++;
                break;
            case 14: // VTK_PYRAMID
                if (nn == 5) {
                    for (int kk = 0; kk < 5; kk++) mm.pyrConn.push_back((int)vv[kk]);
                    mm.pyrRef.push_back(rr);
                } else nUnhandled++;
                break;
            case 13: // VTK_WEDGE (prism)
                if (nn == 6) {
                    for (int kk = 0; kk < 6; kk++) mm.prismConn.push_back((int)vv[kk]);
                    mm.prismRef.push_back(rr);
                } else nUnhandled++;
                break;
            case 12: // VTK_HEXAHEDRON
                if (nn == 8) {
                    for (int kk = 0; kk < 8; kk++) mm.hexConn.push_back((int)vv[kk]);
                    mm.hexRef.push_back(rr);
                } else nUnhandled++;
                break;
            default:
                nUnhandled++;
                break;
        }
        prev = cur;
    }
    if (nUnhandled > 0)
        fprintf(stderr, "VTU: %d cells of unhandled type/arity ignored\n", nUnhandled);

    printf("VTU: %d verts, %d tri, %d quad, %d tet, %d pyr, %d prism, %d hex%s%s%s\n",
           mm.nVerts(), mm.nTri(), mm.nQuad(), mm.nTet(),
           mm.nPyr(), mm.nPrism(), mm.nHex(),
           hasRefs ? " (refs from '" : "",
           hasRefs ? cellRefsName.c_str() : "",
           hasRefs ? "')" : "");
    return 0;
}

static int writeVtu(const char *outputFile, const Mesh &mm)
{
    FILE *fp = fopen(outputFile, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", outputFile); return 1; }

    int64_t nCells = mm.nTri() + mm.nQuad() + mm.nTet() + mm.nPyr() + mm.nPrism() + mm.nHex();

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%lld\">\n",
            mm.nVerts(), (long long)nCells);

    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int ii = 0; ii < mm.nVerts(); ii++)
        fprintf(fp, "          %.17g %.17g %.17g\n", mm.x[ii], mm.y[ii], mm.z[ii]);
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    auto dumpConn = [&](const std::vector<int> &cc, int stride) {
        for (size_t ii = 0; ii < cc.size(); ii += stride) {
            fprintf(fp, "         ");
            for (int kk = 0; kk < stride; kk++) fprintf(fp, " %d", cc[ii + kk]);
            fprintf(fp, "\n");
        }
    };
    dumpConn(mm.triConn,   3);
    dumpConn(mm.quadConn,  4);
    dumpConn(mm.tetConn,   4);
    dumpConn(mm.pyrConn,   5);
    dumpConn(mm.prismConn, 6);
    dumpConn(mm.hexConn,   8);
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    int64_t off = 0;
    auto dumpOff = [&](int stride, int nn) {
        for (int ii = 0; ii < nn; ii++) { off += stride; fprintf(fp, "          %lld\n", (long long)off); }
    };
    dumpOff(3, mm.nTri());
    dumpOff(4, mm.nQuad());
    dumpOff(4, mm.nTet());
    dumpOff(5, mm.nPyr());
    dumpOff(6, mm.nPrism());
    dumpOff(8, mm.nHex());
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    auto dumpType = [&](int tt, int nn) {
        for (int ii = 0; ii < nn; ii++) fprintf(fp, "          %d\n", tt);
    };
    dumpType(5,  mm.nTri());
    dumpType(9,  mm.nQuad());
    dumpType(10, mm.nTet());
    dumpType(14, mm.nPyr());
    dumpType(13, mm.nPrism());
    dumpType(12, mm.nHex());
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    if (nCells > 0) {
        fprintf(fp, "      <CellData Scalars=\"ref\">\n");
        fprintf(fp, "        <DataArray type=\"Int32\" Name=\"ref\" format=\"ascii\">\n");
        auto dumpRef = [&](const std::vector<int> &rr, int nn) {
            for (int ii = 0; ii < nn; ii++)
                fprintf(fp, "          %d\n", rr.empty() ? 0 : rr[ii]);
        };
        dumpRef(mm.triRef,   mm.nTri());
        dumpRef(mm.quadRef,  mm.nQuad());
        dumpRef(mm.tetRef,   mm.nTet());
        dumpRef(mm.pyrRef,   mm.nPyr());
        dumpRef(mm.prismRef, mm.nPrism());
        dumpRef(mm.hexRef,   mm.nHex());
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </CellData>\n");
    }

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    printf("Wrote %s\n", outputFile);
    return 0;
}

// ============================================================
// Dispatch
// ============================================================

enum Format { FMT_UNKNOWN, FMT_STL, FMT_MESHB, FMT_UGRID, FMT_VTU };

static Format detectFormat(const char *path)
{
    if (endsWith(path, ".stl"))                                 return FMT_STL;
    if (endsWith(path, ".meshb") || endsWith(path, ".mesh"))    return FMT_MESHB;
    if (endsWith(path, ".ugrid"))                               return FMT_UGRID;
    if (endsWith(path, ".vtu"))                                 return FMT_VTU;
    return FMT_UNKNOWN;
}

static void printUsage(const char *prog)
{
    printf("Usage: %s input output\n", prog);
    printf("\n");
    printf("Converts between mesh formats. Direction is detected from extensions.\n");
    printf("\n");
    printf("Supported formats:\n");
    printf("  .stl              binary or ASCII STL (surface only)\n");
    printf("  .meshb / .mesh    GMF mesh\n");
    printf("  .ugrid            lb8.ugrid (writes .mapbc companion)\n");
    printf("  .vtu              VTK UnstructuredGrid (ASCII or inline base64 binary,\n");
    printf("                    uncompressed; format=\"appended\" / compressed not supported)\n");
    printf("\n");
    printf("Example: %s model.vtu model.stl\n", prog);
}

int main(int argc, char *argv[])
{
    if (argc < 3) { printUsage(argv[0]); return argc < 2 ? 0 : 1; }

    const char *inFile  = argv[1];
    const char *outFile = argv[2];
    Format inFmt  = detectFormat(inFile);
    Format outFmt = detectFormat(outFile);
    if (inFmt == FMT_UNKNOWN) {
        fprintf(stderr, "Unrecognized input extension: %s\n", inFile);  return 1;
    }
    if (outFmt == FMT_UNKNOWN) {
        fprintf(stderr, "Unrecognized output extension: %s\n", outFile); return 1;
    }

    Mesh mm;
    int rr = 0;
    switch (inFmt) {
        case FMT_STL:   rr = readStl  (inFile, mm); break;
        case FMT_MESHB: rr = readMeshb(inFile, mm); break;
        case FMT_UGRID: rr = readUgrid(inFile, mm); break;
        case FMT_VTU:   rr = readVtu  (inFile, mm); break;
        default: return 1;
    }
    if (rr) return rr;

    fillMissingDerived(mm);

    switch (outFmt) {
        case FMT_STL:   rr = writeStl  (outFile, mm); break;
        case FMT_MESHB: rr = writeMeshb(outFile, mm); break;
        case FMT_UGRID: rr = writeUgrid(outFile, mm); break;
        case FMT_VTU:   rr = writeVtu  (outFile, mm); break;
        default: return 1;
    }
    return rr;
}
