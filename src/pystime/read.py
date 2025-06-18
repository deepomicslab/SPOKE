from __future__ import annotations

import json
import os
import re

from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData, read_text
from h5py import File
from PIL import Image
from scanpy import read_10x_h5, read_10x_mtx, read_csv
from scanpy import logging as logg
from scipy.sparse import csr_matrix

from squidpy._constants._pkg_constants import Key
from squidpy._utils import NDArrayA
from squidpy.datasets._utils import PathLike
from squidpy.read._utils import _load_image, _read_counts

__all__ = ["visium", "vizgen", "nanostring"]


def init_cell_anndata_from_csv(count_fn, cell_meta_fn, gene_meta_fn=None, row_col_order='cell_gene',
                               cell_id_col=None, x_col=None, y_col=None,
                               gene_id_col=None, gene_symbol_col=None,
                               delimiter=','
                              ):
    if row_col_order not in ['gene_cell', 'cell_gene']:
        raise ValueError('Parameter "row_col_order" must be "gene_cell" or "cell_gene", '
                         '"gene_cell" means the row is gene and column is cell, vice versa. '
                         'Currently row_col_order={}.'.format(row_col_order))    
    
    meta_df = pd.read_csv(cell_meta_fn, delimiter=delimiter)

    if 'Unnamed: 0' in meta_df.columns:
        del meta_df['Unnamed: 0']
    
    if cell_id_col not in meta_df.columns:
        raise ValueError('Please specify the column that stands for '
                         'cell or spot id with parameter "cell_id_col". '
                         'Currently cell_id_col={} are not in columns of cell meta file.'.format(cell_id_col))
    meta_df[cell_id_col] = meta_df[cell_id_col].astype(str)

    if x_col not in meta_df.columns:
        raise ValueError('Please specify the column that stands for '
                         'the spatial coordinate x with parameter "x_col". '
                         'Currently x_col={} are not in columns of meta file.'.format(x_col))
        
    if y_col not in meta_df.columns:
        raise ValueError('Please specify the column that stands for '
                         'the spatial coordinate y with parameter "y_col". '
                         'Currently y_col={} are not in columns of meta file.'.format(y_col))
    
    adata = read_csv(count_fn, delimiter=delimiter)
    if row_col_order == 'gene_cell':
        adata = adata.T
    adata.X = csr_matrix(adata.X)
        
    adata.obs_names = adata.obs_names.map(lambda x: x.replace('.','-'))
    print(adata.obs_names)

    if meta_df.shape[0] != adata.shape[0]:
        print('meta', len(meta_df[cell_id_col]))
        print('count', len(adata.obs_names))
        raise ValueError('The number of cell ids in count file and cell meta file not match.')
    
    if not (meta_df[cell_id_col] == adata.obs_names).all():
        print('meta', meta_df[cell_id_col])
        print('count', adata.obs_names)
        raise ValueError('The order of cell ids in count file and cell meta file not match.')
         
    if gene_meta_fn is not None:
        gene_df = pd.read_csv(gene_meta_fn)
        if 'Unnamed: 0' in gene_df.columns:
            del gene_df['Unnamed: 0']
        
        if gene_id_col not in gene_df.columns:
            raise ValueError('Please specify the column that stands for '
                             'gene id with parameter "gene_id_col". '
                             'Currently gene_id_col={} are not in columns of gene meta file.'.format(gene_id_col))        
        
        if not (gene_df[gene_id_col] == adata.var_names).all():
            raise ValueError('The order of gene ids in count file and gene meta file not match.')
        
    adata.X = csr_matrix(adata.X)
    adata.layers['raw'] = adata.X
    
    meta_df.index = meta_df[cell_id_col]
    adata.obs = meta_df
    
    adata.obsm["spatial"] = adata.obs[[x_col, y_col]].to_numpy()
    
    adata.uns['spatial'] = {'library':{}}
    
    if gene_meta_fn is not None:
        if gene_symbol_col is not None:
            gene_df.index = gene_df[gene_symbol_col]
        else:
            gene_df.index = gene_df[gene_id_col]
        adata.var = gene_df
        adata.var['i'] = range(adata.shape[1])
    
    return adata


def visium(
    path: PathLike,
    *,
    counts_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str | None = None,
    load_images: bool = True,
    quality: str | None = 'hires',
    source_image_path: PathLike | None = None,
    **kwargs: Any,
) -> AnnData:
    """
    Read *10x Genomics* Visium formatted dataset.

    In addition to reading the regular *Visium* output, it looks for the *spatial* directory and loads the images,
    spatial coordinates and scale factors.

    .. seealso::

        - `Space Ranger output <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview>`_.
        - :func:`squidpy.pl.spatial_scatter` on how to plot spatial data.

    Parameters
    ----------
    path
        Path to the root directory containing *Visium* files.
    counts_file
        Which file in the passed directory to use as the count file. Typically either *filtered_feature_bc_matrix.h5* or
        *raw_feature_bc_matrix.h5*.
    library_id
        Identifier for the *Visium* library. Useful when concatenating multiple :class:`anndata.AnnData` objects.
    kwargs
        Keyword arguments for :func:`scanpy.read_10x_h5`, :func:`scanpy.read_10x_mtx` or :func:`read_text`.

    Returns
    -------
    Annotated data object with the following keys:
        - :attr:`anndata.AnnData.obs`: observation dataframe with fields 'in_tissue', 'array_row', 'array_col'
        - :attr:`anndata.AnnData.var`: variable dataframe with fields 'gene_ids', 'feature_types'
        - :attr:`anndata.AnnData.obsm` ``['spatial']`` - spatial spot coordinates.
        - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['images']`` - *hires* and *lowres* images.
        - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['scalefactors']`` - scale factors for the spots.
        - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['metadata']`` - various metadata.
    """  # noqa: E501
    path = Path(path)
    adata, library_id = _read_counts(path, count_file=counts_file, library_id=library_id, **kwargs)

    if not load_images:
        return adata

    adata.uns[Key.uns.spatial][library_id][Key.uns.image_key] = {}
    for res in ["fullres", "hires", "lowres"]:
        img_path = path / f"{Key.uns.spatial}/tissue_{res}_image.png"
        if img_path.exists():
            adata.uns[Key.uns.spatial][library_id][Key.uns.image_key][res] = _load_image(img_path) 

    adata.uns[Key.uns.spatial][library_id]["scalefactors"] = json.loads(
        (path / f"{Key.uns.spatial}/scalefactors_json.json").read_bytes()
    )

    tissue_positions_file = (
        path / "spatial/tissue_positions.csv"
        if (path / "spatial/tissue_positions.csv").exists()
        else path / "spatial/tissue_positions_list.csv"
    )

    coords = pd.read_csv(
        tissue_positions_file,
        header=1 if tissue_positions_file.name == "tissue_positions.csv" else None,
        index_col=0,
    )
    coords.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]

    adata.obs = pd.merge(adata.obs, coords, how="left", left_index=True, right_index=True)
    adata.obsm[Key.obsm.spatial] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)

    if source_image_path is not None:
        source_image_path = Path(source_image_path).absolute()
        if not source_image_path.exists():
            logg.warning(f"Path to the high-resolution tissue image `{source_image_path}` does not exist")
        adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(source_image_path)
    
    if quality == "fulres":
        image_coor = adata.obsm["spatial"]
    else:
        scale = adata.uns["spatial"][library_id]["scalefactors"][
            "tissue_" + quality + "_scalef"
        ]
        image_coor = adata.obsm["spatial"] * scale
    
    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]    
    adata.uns["spatial"][library_id]["use_quality"] = quality
    adata.obs['i'] = range(adata.shape[0])
    
    adata.obs["array_row"] = adata.obs["array_row"].astype(int)
    adata.obs["array_col"] = adata.obs["array_col"].astype(int)
    #adata.obsm["spatial"] = adata.obsm["spatial"].astype("int64")
    return adata


def _read_counts(
    path: str | Path,
    count_file: str,
    library_id: Optional[str] = None,
    **kwargs: Any,
) -> Tuple[AnnData, str]:
    path = Path(path)
    if count_file.endswith(".h5"):
        adata: AnnData = read_10x_h5(path / count_file, **kwargs)
        with File(path / count_file, mode="r") as f:
            attrs = dict(f.attrs)
            if library_id is None:
                try:
                    lid = attrs.pop("library_ids")[0]
                    library_id = lid.decode("utf-8") if isinstance(lid, bytes) else str(lid)
                except ValueError:
                    raise KeyError(
                        "Unable to extract library id from attributes. Please specify one explicitly."
                    ) from None

            adata.uns[Key.uns.spatial] = {library_id: {"metadata": {}}}  # can overwrite
            for key in ["chemistry_description", "software_version"]:
                if key not in attrs:
                    continue
                metadata = attrs[key].decode("utf-8") if isinstance(attrs[key], bytes) else attrs[key]
                adata.uns[Key.uns.spatial][library_id]["metadata"][key] = metadata

        return adata, library_id

    if library_id is None:
        raise ValueError("Please explicitly specify library id.")

    if count_file.endswith((".csv", ".txt")):
        adata = read_text(path / count_file, **kwargs)
    elif count_file.endswith(".mtx") or count_file.endswith(".mtx.gz"):
        adata = read_10x_mtx(path / 'filtered_feature_bc_matrix', **kwargs)
    else:
        raise NotImplementedError("TODO")

    adata.uns[Key.uns.spatial] = {library_id: {"metadata": {}}}  # can overwrite
    return adata, library_id


def _load_image(path: PathLike) -> NDArrayA:
    return np.asarray(Image.open(path))