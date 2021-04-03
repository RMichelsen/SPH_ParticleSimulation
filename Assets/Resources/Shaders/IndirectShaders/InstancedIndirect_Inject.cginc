void InjectSetup_float(float3 A, out float3 Out) {
	Out = A;
}

#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
struct IndirectShaderData {
	float4x4 LocalToWorld;
	float4x4 WorldToLocal;
};
#if defined(SHADER_API_GLCORE) || defined(SHADER_API_D3D11) || defined(SHADER_API_GLES3) || defined(SHADER_API_METAL) || defined(SHADER_API_VULKAN) || defined(SHADER_API_PSSL) || defined(SHADER_API_XBOXONE)
uniform StructuredBuffer<IndirectShaderData> _InstancedIndirectShaderData;
#endif

#endif

void SetupProcedural()
{
#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
#ifdef unity_ObjectToWorld
#undef unity_ObjectToWorld
#endif

#ifdef unity_WorldToObject
#undef unity_WorldToObject
#endif
	unity_ObjectToWorld = _InstancedIndirectShaderData[unity_InstanceID].LocalToWorld;
	unity_WorldToObject = _InstancedIndirectShaderData[unity_InstanceID].WorldToLocal;
#endif
}

