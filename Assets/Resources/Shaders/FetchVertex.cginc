StructuredBuffer<float3> vertices;
float visualScale;

void FetchVertexPosition_float(float vertex_id, out float3 Out) {
	Out = vertices[(int)vertex_id * 2] * visualScale;
}
void FetchVertexNormal_float(float vertex_id, out float3 Out) {
	Out = vertices[(int)vertex_id * 2 + 1];
}

