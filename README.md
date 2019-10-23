# sampledBoundedPlane

[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-7-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-7)

## 使用：
### 运行时
- wmake 编译。
- 在 controlDict 文件中使用，参考以下代码设置。
```
functions
{
    surfaceSampling
    {
        type            surfaces;
        functionObjectLibs ( "libZYsampling.so" );
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
                interpolate true;
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
### 运行后处理
- 创建一个文件到 system 文件夹，如 sampledBoundedPlane，内容是：
```
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      sample;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type surfaces;
libs            ("libZYsampling.so");

interpolationScheme cellPoint;

surfaceFormat   vtk;

surfaces
(
    plane1
    {
        type    boundedPlane;
        interpolate true;
        planeType   pointAndNormal;
        pointAndNormalDict
        {
            normal (0 0 1);
            point  (0 0 0) ;
        }
        bounds (0 -0.04 -0.1) (0.08 0.04 0.1);
    }    
);

fields          ( T );


// ************************************************************************* //
```
- 运行 ``postProcess -latestTime -func sampledBoundedPlane``



代码改编自：
https://develop.openfoam.com/Development/OpenFOAM-plus/tree/master/src/sampling/sampledSurface/sampledCuttingPlane
