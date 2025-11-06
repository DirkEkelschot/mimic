#!/usr/bin/env python3

import infinity_plugins as inf
import argparse
import numpy as np
import h5py


def get_clo():
    parser = argparse.ArgumentParser(
        description="Utility to compute the hessian of a scalar field in a snap file and export it as an HDF5 file.")
    parser.add_argument("--mesh", "-m", help="Mesh filename", required=True)
    parser.add_argument("--snap", "-s", help="Snap filename", required=True)
    parser.add_argument(
        "--fieldname", "-f",
        nargs="+",
        default=[],
        help="The name of the field to compute the hessian. Omitting this option computes the Hessian of all fields in the snap file.")
    parser.add_argument(
        "--out", "-o", help="The name of the HDF5 file to output.", default="out.hdf5")
    parser.add_argument(
        "--hessian_plugin", help="T-infinity plugin to use for the Hessian calculator", default="RefinePlugins")
    return parser.parse_args()


def manually_copy_inf_field_to_numpy_array(field: inf.Field):
    out = np.ndarray([field.size(), field.blockSize()])
    for node_id in range(field.size()):
        for eqn, v in enumerate(field.at(node_id)):
            out[node_id, eqn] = v
    return out


@inf.mpi_main
def main():
    mp = inf.getCommunicator()
    opts = get_clo()
    print(f"loading mesh: {opts.mesh}")
    mesh = inf.loadMesh(mp, opts.mesh)
    snap = inf.Snap(opts.snap, mesh, mp)

    gradient_calculator = inf.GradientCalculator(
        inf.get_plugin_library_path(), opts.hessian_plugin, mp, mesh, "{}")

    available_fields = snap.listFields()
    if opts.fieldname:
        fields_to_compute = opts.fieldname
    else:
        fields_to_compute = available_fields

    for f in fields_to_compute:
        assert (f in available_fields)

    with h5py.File(opts.out, "w") as f:
        for name in fields_to_compute:
            print(f"Computing hessian of scalar field <{name}>")
            field = snap.field(name)
            f[name] = field.to_array()
            hessian = gradient_calculator.calcHessian(field)
            f[f"{name}-hessian"] = manually_copy_inf_field_to_numpy_array(inf.reflectHessian(hessian))
        print(f"writing file <{opts.out}> to disk... ", end="")
    print("Done")


if __name__ == "__main__":
    main()
 
