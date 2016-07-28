function matconvnet_setup_GPU(varargin)

deepagpaths();
vl_setupnn();

opts.useGpu = true;
opts.enableCudnn = true;
opts.cudnnRoot = '/meleze/gpudata/SHARED/cudnn/cudnnR4-prod/';
opts.verbose = false;
opts = vl_argparse(opts, varargin) ;

try
  vl_nnconv(single(1),single(1),[]) ;
catch
  warning('VL_NNCONV() does not seem to be compiled. Trying to compile it now.') ;
  vl_compilenn('enableGpu', opts.useGpu, 'verbose', opts.verbose) ;
end

if opts.useGpu
  try
    vl_nnconv(gpuArray(single(1)),gpuArray(single(1)),[]) ;
  catch
    vl_compilenn('enableGpu', opts.useGPU, 'enableCudnn', opts.enableCudnn, 'cudnnRoot', opts.cudnnRoot, 'verbose', opts.verbose);
    warning('GPU support does not seem to be compiled in MatConvNet. Trying to compile it now') ;
  end
end
