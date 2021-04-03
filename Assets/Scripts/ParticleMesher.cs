using System.Text;
using System.IO;
using Unity.Mathematics;
using UnityEngine;
using System.Collections.Generic;

public static class ParticleMesher {
    public static void WriteParticlesToDisk(float3x2[] vertices, string path) {

        StringBuilder vertexString = new StringBuilder();
        StringBuilder normalString = new StringBuilder();
        StringBuilder faceString = new StringBuilder();
        faceString.Append("f ");

        for (int i = 0; i < vertices.Length; ++i) {
            float3x2 vertex = vertices[i];
            float3 P = vertex.c0;
            float3 N = math.normalize(vertex.c1);
            vertexString.AppendFormat("v {0} {1} {2}\n", P.x, P.y, P.z);
            normalString.AppendFormat("vn {0} {1} {2}\n", N.x, N.y, N.z);
            faceString.AppendFormat("{0}//{1} ", i + 1, i + 1);

            if((i + 1) % 3 == 0) {
                if(i == vertices.Length - 1) {
                    faceString.Append("\n");
                }
                else {
                    faceString.Append("\nf ");
                }
            }
        }

        using (StreamWriter sw = File.CreateText(path)) {
            sw.Write(vertexString);
            sw.Write(normalString);
            sw.Write(faceString);
        }
    }
}
