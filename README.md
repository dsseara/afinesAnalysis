# AFiNEs Analysis tools

A suite of matlab functions and scripts for analyzing output from AFiNEs simulations.

In general, the workflow looks something like the following:
First, follow the instructions in the first link below to download and run AFiNEs on the Grace cluster. While that's going on, clone this repository into a location of your choosing:
```bash
git clone https://github.com/dsseara/afinesAnalysis.git
```
*Make sure that matlab can find this folder by adding it to your path in your startup.m file*
Once you have the output in the appropriate data/ and txt_stack/ folders, cd to their parent folder within matlab and enter 
```MATLAB
run('read_data.m');
run('velInterp.m');
run('divVelocity.m');
run('divStats.m');
```
`read_data.m` will produce a file called `simdata.mat` that contains positions over time of the actin, motors, and crosslinkers in `adata`, `mdata`, and `pdata` respectively (the 'p' in `pdata` refers to the fact that AFiNEs calls crosslinkers 'passive motors'). `adata` is an Nx4xM matrix, where N is the number of actin beads and M is the number of time steps. The second dimension is of the form  [x_position y_position link_length filament_id]. There is also a `params` structure created. For example, for a simulation where 1000 frames are taken over 100 seconds in a domain of size 100um x 100um
```MATLAB
>> params

params = 

     timestep: [1x1000 double]
         anum: [1x1000 double]
       xRange: [-50 50]
       yRange: [-50 50]
            L: 100
           dt: 0.1000
        npoly: [1000x1 double]
     numBeads: 10999
          dim: 4
    numFrames: 1000
         mnum: [1x1000 double]
         pnum: [1x1000 double]
```


## AFiNEs

Find the [software][1] with instructions, and the [paper][2] for more details about how it's written.


[1]: https://github.com/shilobanerjee/AFiNeS
[2]: https://arxiv.org/abs/1609.05202
