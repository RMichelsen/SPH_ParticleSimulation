#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
StructuredBuffer<float3> colors;

void GetColor_float(out float3 Out) {
	Out = colors[unity_InstanceID];
}
#else
void GetColor_float(out float3 Out) {
	Out = float3(0.0f, 0.0f, 0.0f);
}
#endif
