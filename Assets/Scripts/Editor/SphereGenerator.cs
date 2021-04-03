using UnityEngine;
using UnityEditor;
using System.Collections.Generic;
using System;
using UnityEngine.Rendering;

public class SphereGenerator : MonoBehaviour
{
    private const float S = 0.5257311121191336f;
    private const float T = 0.85065080835204f;
    private const float U = 0.0f;

    private static int DivideVertex(List<Vector3> vertices, Dictionary<Tuple<int, int>, int> edgeLookup, int v0, int v1) {
        Tuple<int, int> pair = new Tuple<int, int>(v0, v1);
        if(v0 > v1) {
            pair = new Tuple<int, int>(v1, v0);
        }
        int element = 0;
        if(!edgeLookup.ContainsKey(pair)) {
            element = (int)vertices.Count;
            edgeLookup.Add(pair, (int)vertices.Count);
            vertices.Add(Vector3.Normalize(vertices[v0] + vertices[v1]));
        }
        else {
            element = edgeLookup[pair];
        }

        return element;
    }

    private static List<int> Subdivide(List<Vector3> vertices, List<int> indices) {
        Dictionary<Tuple<int, int>, int> edgeLookup = new Dictionary<Tuple<int, int>, int>();
        List<int> resultingIndices = new List<int>();

        int initialSize = indices.Count;
        for(int i = 0; i < initialSize; i += 3) {
            int mid1 = DivideVertex(vertices, edgeLookup, indices[i + 0], indices[i + 1]);
            int mid2 = DivideVertex(vertices, edgeLookup, indices[i + 1], indices[i + 2]);
            int mid3 = DivideVertex(vertices, edgeLookup, indices[i + 2], indices[i + 0]);

            resultingIndices.AddRange(new List<int>() { indices[i + 0], mid1, mid3 });
            resultingIndices.AddRange(new List<int>() { indices[i + 1], mid2, mid1 });
            resultingIndices.AddRange(new List<int>() { indices[i + 2], mid3, mid2 });
            resultingIndices.AddRange(new List<int>() { mid1, mid2, mid3 });
        }

        return resultingIndices;
    }

    [MenuItem("Mesh Generation/Create Sphere")]
    public static void CreateSphere()
    {
        List<Vector3> vertices = new List<Vector3>() {
            new Vector3(-S,  U,  T), new Vector3(S,  U,  T),  new Vector3(-S,  U, -T), new Vector3( S,  U, -T),
            new Vector3(U,  T,  S),  new Vector3(U,  T, -S),  new Vector3( U, -T,  S), new Vector3( U, -T, -S),
            new Vector3(T,  S,  U),  new Vector3(-T,  S,  U), new Vector3( T, -S,  U), new Vector3(-T, -S, U)
        };
        List<int> indices = new List<int> {
            0, 4, 1,  0, 9, 4,  9, 5, 4,  4, 5, 8,  4, 8, 1,
            8, 10, 1, 8, 3, 10, 5, 3, 8,  5, 2, 3,  2, 7, 3,
            7, 10, 3, 7, 6, 10, 7, 11, 6, 11, 0, 6, 0, 1, 6,
            6, 1, 10, 9, 0, 11, 9, 11, 2, 9, 2, 5,  7, 2, 11
        };

        for(int i = 0; i < 2; i++) {
            indices = Subdivide(vertices, indices);
        }

        Mesh mesh = new Mesh();
        mesh.indexFormat = IndexFormat.UInt16;
        mesh.vertices = vertices.ToArray();
        mesh.triangles = indices.ToArray();
        mesh.RecalculateNormals();

        AssetDatabase.CreateAsset(mesh, "Assets/Resources/Meshes/Sphere.asset");
        AssetDatabase.SaveAssets();
    }
}