# sampledBoundedPlane

## 使用：
- wmake 编译。
- 在 controlDict 文件中加入 libs ("libZYsampling.so");
- 在 controlDict 文件中使用，参考以下代码设置。
```
functions
{
    surfaceSampling
    {
        type            surfaces;
        functionObjectLibs ( "libsampling.so" );
        enabled         true;
        writeControl    timeStep;
        writeInterval   1;
        interpolationScheme cellPointFace;
        surfaceFormat   vtk;
        fields          (T);
        surfaces
        (
            plane1
            {
                type    boundedPlane;
                planeType   pointAndNormal;
                pointAndNormalDict
                {
                    normal (0 1 0);
                    point  (0 0 0) ;
                }
                bounds (0 -0.006 -0.1) (0.15 0.006 0.1);
            }    
        );
    }
}
```
