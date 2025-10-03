# trigeminal

This repo allows you to run a pipeline to track the trigeminal nerves.
There are two scripts that need to be run one after the other.

You will need scilpy 2.2.0 and ANTs to be installed in order to run these scripts.

- trigeminal_first_order.sh
- trigeminal_second_order.sh

## Input structure

```
    [input]
    ├── sub-01
    │   ├── freesurfer
    │   │   └─── aparc.DKTatlas+aseg.mgz
    │   ├── sub-01__fa.nii.gz
    │   ├── sub-01__fodf.nii.gz
    │   └── sub-01__t1_warped.nii.gz
    │
    ├── S2
    .
    .
```

## Examples

### First order
```
trigeminal_first_order.sh \
    -s path/to/input \
    -m ROIs_clean \
    -o my_results \
    -t 10 \ # Number of threads
    -g true # Use GPU
```

### Second order
```
trigeminal_second_order.sh \
    -s path/to/input \
    -m ROIs_clean \
    -o my_results \
    -t 10 \ # Number of threads
    -g true # Use GPU
```
