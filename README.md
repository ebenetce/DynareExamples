<a name="T_5442E3AB"></a>
# <span style="color:rgb(213,80,0)">Dynare examples</span>
<a name="beginToc"></a>
## Table of Contents
&emsp;[ramst.mod](#H_2765B1D4)
 
<a name="endToc"></a>
<a name="H_2765B1D4"></a>
## ramst.mod

An elementary real business cycle (RBC) model, simulated in a deterministic setup.

<pre>
cd ./models/ramst
dynare ramst.mod
</pre>
<a name="H_924133BF"></a>
## <samp>example1.mod</samp> [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ebenetce/DynareExamples&file=DynareExamples.mlx&line=1)

An example of a small RBC model in a stochastic setup, presented in *Collard (2001)* (see the file <samp>guide.pdf</samp> which comes with Dynare).

<pre>
cd ./models/example1
dynare example1.mod
</pre>
<a name="H_057F809F"></a>
## <samp>example2.mod</samp> [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ebenetce/DynareExamples&file=DynareExamples.mlx&line=2)

An example of a small RBC model in a stochastic setup, presented in *Collard (2001)* (see the file <samp>guide.pdf</samp> which comes with Dynare).

<pre>
cd ./models/example2
dynare example2.mod
</pre>
<a name="H_32D5D587"></a>
## <samp>example3.mod</samp> [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ebenetce/DynareExamples&file=DynareExamples.mlx&line=3)

A small RBC model in a stochastic setup, presented in *Collard (2001)*. The steady state is solved analytically using the <samp>steady_state_model</samp> block (see [**<samp>steady_state_model</samp>**](https://www.dynare.org/manual/the-model-file.html#steady_state_model)).

<pre>
cd ./models/example3
dynare example3.mod
</pre>
<a name="H_1DA975A5"></a>
## <samp>fs2000.mod</samp>

A cash in advance model, estimated by *Schorfheide (2000)*. The file shows how to use Dynare for estimation.

<pre>
cd models\fs2000
dynare fs2000.mod
</pre>
<a name="H_FBB805EB"></a>
## <samp>fs2000_nonstationary.mod</samp>

The same model than <samp>fs2000.mod</samp>, but written in non-stationary form. Detrending of the equations is done by Dynare.

<pre>
cd models\fs2000_nonstationary
dynare fs2000_nonstationary.mod
</pre>
<a name="H_C8E8493C"></a>
## <samp>bkk.mod</samp>

Multi-country RBC model with time to build, presented in *Backus, Kehoe and Kydland (1992)*. The file shows how to use Dynare’s macro processor.

<pre>
cd ./models/bkk
dynare bkk.mod
</pre>
<a name="H_0D47BFA4"></a>
## <samp>agtrend.mod</samp>

Small open economy RBC model with shocks to the growth trend, presented in *Aguiar and Gopinath (2004)*.

<pre>
cd ./models/agtrend
dynare agtrend.mod
</pre>
<a name="H_BA98565B"></a>
## <samp>Gali_2015.mod</samp>

Basic New Keynesian model of *Galí (2015)*, Chapter 3 showing how to i) use “system prior”-type prior restrictions as in *Andrle and Plašil (2018)* and ii) run prior/posterior-functions.

<pre>
cd ./models/Gali_2015/
dynare Gali_2015.mod
</pre>
<a name="H_D48BD1C1"></a>
## <samp>NK_baseline.mod</samp>

Baseline New Keynesian Model estimated in *Fernández-Villaverde (2010)*. It demonstrates how to use an explicit steady state file to update parameters and call a numerical solver.

<pre>
cd ./models/NK_baseline
dynare NK_baseline.mod
</pre>
<a name="H_6AC90C0E"></a>
## <samp>Occbin_example.mod</samp>

RBC model with two occasionally binding constraints. Demonstrates how to set up Occbin.

<pre>
cd ./models/Occbin
dynare Occbin_example.mod
</pre>
<a name="H_CF1D1EEA"></a>
## <samp>Ramsey_Example.mod</samp>

File demonstrating how to conduct optimal policy experiments in a simple New Keynesian model either under commitment (Ramsey) or using optimal simple rules (OSR)

<pre>
cd ./models/Ramsey_Example/
dynare Ramsey_example.mod
</pre>
<a name="H_2198765D"></a>
## <samp>Ramsey_steady_file.mod</samp>

File demonstrating how to conduct optimal policy experiments in a simple New Keynesian model under commitment (Ramsey) with a user-defined conditional steady state file

<pre>
cd ./models/Ramsey_steady_file/
dynare Ramsey_steady_file.mod
</pre>
