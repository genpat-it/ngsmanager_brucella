# NGSMANAGER

`ngsmanager` is a bioinformatics best-practice analysis pipeline built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Prerequisites

At this stage it is possibile to run a workflow (pipeline or module) logging into `tst-nextflow`, a server where:

- nextflow and its dependencies are configured and installed
- a git client is installed
- it is possibile to access this git repo
- mount-points for external resources (like databases) are configured
- it is possibile to access the internal docker repo to get all the required docker images

:warning: make sure you have `gtc-devsrv:3000` configured in your ~/.nextflow/scm as gitea provider.
If not, please add to that file:
```
providers {
    gtc-devsrv {
        server = 'http://gtc-devsrv:3000'
        endpoint = 'http://gtc-devsrv:3000/api/v1'
        platform = 'gitea'
    }
}
```

## How to run a pipeline or module

The typical command for running a pipeline or a module is as follows:

```
nextflow run -latest http://gtc-devsrv:3000/bioinfo/ngsmanager -entry <pipeline or module entrypoint> --cmp <sample identifier> --dt <date of exam, format: 'DTyyMMdd'>  --ds <sample distribution code>
```
where:
- `cmp` is the identifier of the sample (e.g. '2020.TE.90974.1.5')
- `dt` is the date of the exam whose output will be considered as the input for the step/module/pipeline specified. (e.g. 'DT200310')
- `ds` is the sample distribution code  (e.g. 'DS10570914')

Therefore, to execute the trimming for the raw reads related to sample '2020.TE.90974.1.5' with DS 'DS10570914' as examined in date 'DT200310' we could run:

```
nextflow run -latest http://gtc-devsrv:3000/bioinfo/ngsmanager --cmp 2020.TE.90974.1.5 --dt DT200310 --ds DS10570914 -entry s_1PP_trimming
```

Note that the pipeline will create the following files in your current directory:

```
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
```
### Common Inputs Parameters 

- `--outdir`: directory where the results will be saved (default: `./results`)
- `--inputDir`: directory used to get the sample data (default: `/bioinfonas/cmp`)
- `--tracedir`: directory used to get the execution info (default: ` "${params.outdir}/pipeline_info"`)

### Single sample execution - Parameters 

- `--cmp`: identifier of the sample (e.g. '2020.TE.90974.1.5')
- `--dt`: date of the exam whose output will be considered as the input for the step/module/pipeline specified. (e.g. 'DT200310')
- `--ds`: the distribution code of the sample (e.g. 'DS10570914')

### Multi sample execution (from SampleSheet) - Parameters 

- `--runPath`: path where sequencer results are stored
- `--sampleSheet`: SampleSheet CSV file, containing metadata of samples to be analyzed (cmp, ds, dt and much more)