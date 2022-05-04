# STOne-Foveation
Post-Acquisition Single Pixel Compressive Foveation using the STOne Transform

![Fig_1_snip](https://user-images.githubusercontent.com/44003292/166619966-c6111df1-d469-4880-a91b-38db367d1819.PNG)

In standard single-pixel foveation methods, the foveated regions need to be determined before any measurements are acquired. This method, however, allows for the determination of the foveated region after all measurements have been acquired, is entirely non-adaptive, and requires only that random STOne patterns be used to sample the scene.

This is accomplished by taking advantage of the downsampling property of the STOne matrix. If a high-resolution STOne pattern is downsampled, then a low-resolution STOne pattern is obtained. A corollary to this is that any measurement acquired with a high-resolution STOne pattern is approximately equal to a measurement acquired with the corresponding low-resolution pattern. Then by sampling the scene with high-resolution STOne patterns but reconstructing using foveated STOne patterns, foveated reconstructions can be easily achieved. Because foveated reconstructions recover fewer overall pixels, an increase in reconstrution time is seen.

![Fig_2_snip](https://user-images.githubusercontent.com/44003292/166619927-f9988804-1692-4e98-a494-4db76988ec9b.PNG)

To obtain the hyperspectral video data, obtain the "Ball" test data from here: https://bearshng.github.io/mht/.

Hyperspectral static data available upon request.

If you use this code, please cite the original paper:

Giljum, A.T. and Kelly, K.F., "On-the-Fly Compressive Single-Pixel Foveation using the STOne Transform," _Opt. Expr._, 2022.
