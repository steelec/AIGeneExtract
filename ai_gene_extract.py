


def read_donor_data(data_dir,filter_noise=True):
    ## Adpated from code that K. Gorgolewski wrote (alleninf) TODO: get full link to github repo
    # expanded to include well data
    # reads allen database donor data (whatever is in the dir)
    import pandas as pd
    from glob import glob
    import os

    donor_ids = [path.split(os.path.sep)[-2] for path in glob(os.path.join(data_dir, "*", "MicroarrayExpression.csv"))]
    print "Data directory contains the following donors: %s" % ", ".join(donor_ids)
    main_df = "empty"
    for donor_id in donor_ids:
        print "Reading data from donor %s"%donor_id
        sample_locations = pd.read_csv(os.path.join(data_dir, donor_id, 'SampleAnnot.csv'))
        df = pd.DataFrame([donor_id] * sample_locations.shape[0], columns=["donor_id"]) #column of donor_id
        df_well_ids = sample_locations[['well_id',"slab_type"]] #keep slab type so that we can filter on it
        df = pd.concat([df,df_well_ids],axis=1,ignore_index=False)
        expression_data = pd.read_csv(os.path.join(data_dir, donor_id, 'MicroarrayExpression.csv'), header=None, index_col=0)
        if filter_noise: #filter to only include data above the noise floor
            good_data = pd.read_csv(os.path.join(data_dir, donor_id, 'PACall.csv'), header=None, index_col=0)
            expression_data[good_data == 0] = np.nan
        expression_data.columns = range(expression_data.shape[1])
        df = pd.concat([df, expression_data.T], axis=1, ignore_index=False)
        if isinstance(main_df, str):
            main_df = df
        else:
            main_df = pd.concat([main_df, df], ignore_index=True)
    return main_df #now formatted by row,col sample_well_id,probe_id


def mm2vox(aff,pts):
    import nibabel as nb
    import numpy as np
    #convert xyz coords from mm to voxel space coords
    return (nb.affines.apply_affine(np.linalg.inv(aff),pts)).astype(int)
def vox2mm(aff,pts):
    import nibabel as nb
    #convert from voxel coords space back to mm space xyz
    return nb.affines.apply_affine(aff,pts)

def get_probe_ids(gene_symbols, df_probes2genes_mapping,
                  gene_sybmol_col_name='gene_symbol',
                  probe_col_name='probe_id'):
    # extract the probe ids when provided with single/multiple gene_symbols and the probe to genes mapping

    df_probes2genes_mapping = df_probes2genes_mapping.reset_index()
    df_probes2genes_mapping = df_probes2genes_mapping.set_index(gene_sybmol_col_name)

    if isinstance(gene_symbols, (
    int, long)):  # cast to list if we only provided the single gene, since casts to int if <2 probes are returned
        gene_symbols = [gene_symbols]
    elif isinstance(gene_symbols, (str)):
        gene_symbols = [gene_symbols]

    probe_ids = df_probes2genes_mapping.loc[gene_symbols][probe_col_name]
    probe_ids = probe_ids.reset_index()
    return probe_ids


def get_gene_expression_old(probe_id, df_donor_data,
                            donor_col_name='donor_id', well_col_name='well_id'):
    # pull gene expression across wells for the gene (probe_id) provided
    # returns dataframe with donor_id, well_id, and expression values for this probe
    # only does a single probe_id, so don't try more :-(

    # donor_data has columns for each well_id and for each gene probe number
    # probe_ids are numbers, so cast to string to select the column
    #    if isinstance(probe_id,(int,long)):
    #        probe_id = str(probe_id)
    #        print("yo dude?")
    #    else:
    #        probe_id = [str(x) for x in probe_id]
    gene_expression = df_donor_data[[donor_col_name, well_col_name, probe_id]]
    # cols = gene_expression.columns
    # gene_expression.columns = cols[0:2]+'gene_expression'
    gene_expression.rename(columns={probe_id: 'expression'}, inplace=True)
    return gene_expression


def get_gene_expression(probe_ids, df_donor_data,
                        donor_col_name='donor_id', well_col_name='well_id',
                        gene_symbol_name="gene_symbol", mean_data=True):
    # pull gene expression across wells for the gene (probe_ids) provided
    # returns dataframe with donor_id, well_id, and average expression values for these probes per well

    if isinstance(probe_ids, (int, long)):
        gene_expression = df_donor_data[[donor_col_name, well_col_name, probe_ids]].copy()
        gene_expression.rename(columns={probe_id: 'expression'}, inplace=True)
    else:
        print("Averaging {0} expression across {1} probes".format("oops", len(list(probe_ids))))
        cols = list([donor_col_name])
        cols.append(well_col_name)
        cols.extend(list(probe_ids))
        gene_expression = df_donor_data[cols].copy()
        gene_expression = gene_expression.assign(expression=gene_expression[list(probe_ids)].mean(axis=1))
    return gene_expression


def get_MNI_coords(gene_expression, well2MNI_mapping, well_col_name='well_id'):
    import pandas as pd
    well2MNI_mapping = well2MNI_mapping.reset_index()
    well2MNI_mapping = well2MNI_mapping.set_index(well_col_name)
    MNI_coords = well2MNI_mapping.loc[gene_expression[well_col_name]].reset_index()  # list of coords
    # ordering is preserved, so can concatenate (merging duplicates rows due to repeated well_ids if >1 gene selected)
    # this works, but creates a duplicate
    # df.T.drop_duplicates().T
    return pd.concat([gene_expression, MNI_coords], axis=1)  # this ignores any ordering, assumes all ordered correctly

def plot_gene_expression(gene_expression_coords, MNI_template=None,
                         well_col_name='well_id', expression_col_name='expression',
                         gene_symbol=None, gene_symbol_col_name = "gene_symbol",
                         flip_coords_to_left_hemisphere=False,
                         add_vox_xyz=None, smoothing_kernel=5,
                         zscore_plotting_data=False, black_bg=False,
                         nan_val=0, plot_smoothed_data=False, donor_id=None,
                         out_dir=None):
    # plot data with nilearn glass brain, weeeeeeee
    # flip_coords_to_left_hemisphere to put all datapoints into single hemi in case data is too sparse
    # TODO:  average across all genes when gene symbol set to None

    import nibabel as nb
    import numpy as np
    from nilearn import plotting
    from nilearn.image import smooth_img
    from nilearn.image import math_img
    import os

    if gene_symbol is not None:
        gene_expression_coords = gene_expression_coords[gene_expression_coords[gene_symbol_col_name] == gene_symbol]
    else:
        gene_symbol = "XXX"
        gene_expression_coords = gene_expression_coords.T.drop_duplicates().T #remove duplicated column
        print("This is not currently functional. Results are not correct.")
        #TODO: calculate the average if you feel like it
    # from scipy.stats import zscore

    if MNI_template is None: #assume that we can find it here :-/
        MNI_template = "/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz"

    # get MNI template
    img = nb.load(MNI_template)
    MNI_d = img.get_data()
    MNI_aff = img.affine
    h = img.header

    df_orig = gene_expression_coords.copy()
    genes = df_orig[gene_symbol_col_name].unique()
    # XXX loop over each of the genes
    # to be able to handle that case
    # map voxel coordinates
    gene_loc_d = np.zeros_like(MNI_d).astype(np.float32)  # will be filled with points from our wells!

    if out_dir is None:
        out_dir = os.path.curdir

    donor_tag = "all"
    if donor_id is not None:
        try:
            gene_expression_coords = gene_expression_coords[gene_expression_coords['donor_id'] == donor_id]
            donor_tag = donor_id.replace(".", "p")
        except:
            print(
            "Tried to filter by donor_id == {0} but failed\n --> please check if this donor_id is included in the dataset".format(
                donor_id))
            print(list(gene_expression_coords['donor_id'].unique()))
            return list(gene_expression_coords['donor_id'].unique())
    if flip_coords_to_left_hemisphere:
        print(
        "You chose to flip all data to the left hemisphere, there is NO guarantee that there will not\nbe overlap that is misrepresented in the figure.")
        gene_expression_coords.iloc[:, -3] = np.abs(
            gene_expression_coords.iloc[:, -3]) * -1  # XXX ALWAYS MUST BE LAST THREE

    if zscore_plotting_data:
        gene_expression_coords[expression_col_name] = (gene_expression_coords[expression_col_name].values - np.nanmean(
            gene_expression_coords[expression_col_name].values)) / np.nanstd(
            gene_expression_coords[expression_col_name].values)
        # print( zscore(gene_expression_coords[expression_col_name].values))

    print("number of wells in data: {0}".format(gene_expression_coords.shape[0]))
    for index, row in gene_expression_coords.iterrows():
        coord = row.values[-3:]  # XXX not the best way to do this :-/
        coord = mm2vox(MNI_aff, coord)  # convert to voxel space for plotting
        if np.isnan(row[expression_col_name]):
            val = nan_val
        else:
            val = row[expression_col_name]
        gene_loc_d[coord[0], coord[1], coord[2]] = val
        # fill it out a bit for display XXX temporary XXX, need to implement something similar to Chris'
        # or keep as is and work via closeness to volumetric structure or surface
        if (add_vox_xyz is not None) and (add_vox_xyz > 0):
            gene_loc_d[coord[0] - add_vox_xyz:coord[0] + add_vox_xyz + 1,
            coord[1] - add_vox_xyz:coord[1] + add_vox_xyz + 1, coord[2] - add_vox_xyz:coord[2] + add_vox_xyz + 1] = val

    out_img = nb.Nifti1Image(gene_loc_d, MNI_aff, header=h)
    out_img.set_data_dtype('float32')
    out_img.header['cal_min'] = np.min(gene_loc_d)
    out_img.header['cal_max'] = np.max(gene_loc_d)
    out_img.update_header()
    img_fname = os.path.join(out_dir, "allen_well_locations_expression_" + donor_tag + "_" + gene_symbol + ".nii.gz")
    print(img_fname)
    nb.save(out_img, img_fname)

    if plot_smoothed_data:
        out_img_s = smooth_img(out_img, smoothing_kernel)
        # out_img_s = math_img("img * 100",img = out_img_s)
        out_img.set_data_dtype('float32')
        out_img_s.header['cal_min'] = np.min(out_img_s.get_data())
        out_img_s.header['cal_max'] = np.max(out_img_s.get_data())
        out_img_s.update_header()
        img_fname = os.path.join(out_dir,
                                 "allen_well_locations_expression_" + donor_tag + "_" + gene_symbol + "_s" + str(
                                     smoothing_kernel) + ".nii.gz")
        print(img_fname)
        nb.save(out_img_s, img_fname)
    print("Mapping of wells to voxels saved to nii.gz file.")
    if np.min(gene_loc_d) < 0:
        plotting.plot_glass_brain(out_img, black_bg=black_bg, display_mode='lyrz',
                                  colorbar=True, title="Allen well locations zscore(" + " [" + gene_symbol + "])",
                                  plot_abs=False)
        if plot_smoothed_data:
            plotting.plot_glass_brain(out_img_s, black_bg=black_bg, display_mode='lyrz',
                                      colorbar=True, title="Allen well locations (zscore smoothed " + str(
                    smoothing_kernel) + "mm fwhm)" + " [" + gene_symbol + "]",
                                      plot_abs=False)
    else:
        plotting.plot_glass_brain(out_img, black_bg=black_bg, display_mode='lyrz',
                                  colorbar=True, title="Allen well locations (expression)" + " [" + gene_symbol + "]")
        if plot_smoothed_data:
            plotting.plot_glass_brain(out_img_s, black_bg=black_bg, display_mode='lyrz',
                                      colorbar=True, title="Allen well locations (smoothed " + str(
                    smoothing_kernel) + "mm fwhm)" + " [" + gene_symbol + "]")


def get_gene_expression_summary_in_mask(nii_mask, gene_expression_coords, nan_val=0,
                                        expression_col_name="expression", gene_symbol_col_name = "gene_symbol",
                                        summary_type="average", donor_col_name='donor_id', verbose=False):
    # compute summary metrics in provided mask, computes the summary per donor_id
    # mask is 3d, but can contiain multiple ids, I think I skip zero?
    # currently assumes that ALL DATA IS FOR THE SAME GENE!!!
    # TODO: iterate over multiple genes as well
    # TODO: zfill, take index mapping to names

    import nibabel as nb
    import numpy as np
    import pandas as pd

    img = nb.load(nii_mask)
    mask_d = img.get_data()
    MNI_aff = img.affine
    vals = np.unique(mask_d)
    mask_ids = np.sort(vals[np.where(vals)]).astype(int)  # nonzero index values, sorted

    if len(mask_ids) == 1:
        mask_ids = [mask_ids]
    # convert gene expression data to voxel space for all locations for each donor, then calculate the summary
    donor_ids = gene_expression_coords[donor_col_name].unique()

    # create a dataframe to store the results
    cols = [donor_col_name, gene_symbol_col_name]
    gene_symbols = np.unique(gene_expression_coords[gene_symbol_col_name])
    gene_symbol = gene_symbols[0] #TODO: iterate over multiple genes, CURRENTLY ONLY TAKES A SINGLE ONE

    for mask_id in mask_ids:
        cols.append("mask_id_" + str(mask_id))
    index = np.arange(0, len(donor_ids))
    df_res = pd.DataFrame(columns=cols, index=index)
    df_res[donor_col_name] = donor_ids
    df_res = df_res.set_index(donor_col_name)
    for donor_id in donor_ids:
        gene_loc_d = np.zeros_like(mask_d)

        # put the data into the same 3d form that the mask is in
        print("Creating the 3d expression data from donor {0}".format(donor_id))
        for index, row in gene_expression_coords[gene_expression_coords[donor_col_name] == donor_id].iterrows():
            coord = row.values[-3:]  # XXX not the best way to do this :-/
            coord = mm2vox(MNI_aff, coord)  # convert to voxel space for plotting
            if np.isnan(row[expression_col_name]):
                expression_val = nan_val
            else:
                expression_val = row[expression_col_name]
                gene_loc_d[coord[0], coord[1], coord[2]] = expression_val
                #   print(expression_val)
        if verbose:
            print("  Extracting mask data")
            print("  mask_id (mean): "),
        for mask_id in mask_ids:
            if verbose:
                print(mask_id),
            single_mask = np.logical_and(np.logical_not(np.isnan(mask_d)), mask_d == mask_id)
            region = np.logical_and(gene_loc_d, single_mask)  # all voxels that have values and are within the mask
            if summary_type is "average":
                expression_summary_val = gene_loc_d[region].mean()
            if verbose:
                print("({:.2f})".format(expression_summary_val)),
            # return df_res
            df_res.loc[donor_id, "mask_id_" + str(mask_id)] = expression_summary_val
            df_res.loc[donor_id, gene_symbol_col_name] = gene_symbol
        if verbose:
            print("")
    return df_res.reset_index()


def get_gene_expression_multi(df_probe_ids, df_donor_data,
                              donor_col_name='donor_id', well_col_name='well_id',
                              probe_col_name='probe_id', gene_symbol_col_name="gene_symbol",
                              mean_data=True):
    # pull gene expression across wells for the gene (df_probe_ids) provided
    # returns dataframe with donor_id, well_id, and average expression values for these probes per well
    # does this for multiple genes, returns all expression values and single column with mean expression value

    import pandas as pd
    all_cols = list([donor_col_name])
    all_cols.append(gene_symbol_col_name)
    all_cols.append(well_col_name)
    if not isinstance(df_probe_ids[probe_col_name], (int, long)):
        all_cols.extend(list(df_probe_ids[probe_col_name]))
    else:
        all_cols.append(df_probe_ids[probe_col_name])
    gene_expression = pd.DataFrame(columns=all_cols)
    gene_symbols = df_probe_ids[gene_symbol_col_name]
    if isinstance(gene_symbols, str):
        gene_symbols = [gene_symbols]
    else:
        gene_symbols = gene_symbols.unique()
    for gene_symbol in gene_symbols:
        print("Collecting gene expression: [{0}]".format(gene_symbol))
        df_gene = df_probe_ids[df_probe_ids[gene_symbol_col_name] == gene_symbol]
        if isinstance(df_gene[probe_col_name], (int, long)):  # XXX this case hasn't been tested as of yet
            gene_expression = df_donor_data[[donor_col_name, well_col_name, df_gene[probe_col_name]]].copy()
            gene_expression.insert(1, gene_symbol_col_name, gene_symbol)
            gene_expression.insert(-1, "expression", df_gene[probe_col_name])  # just fill with the value of the probe
        # gene_expression.rename(columns={probe_id:'expression'+"_"+gene_symbol},inplace=True) #this is the original data, just renamed since there is no average
        else:
            print("  Averaging expression across {0} probes".format(len(list(df_gene[probe_col_name]))))
            cols = list([donor_col_name])
            cols.append(well_col_name)
            cols.extend(list(df_gene[probe_col_name]))
            df_donor_data_gene_probes = df_donor_data[cols].copy()
            df_donor_data_gene_probes.insert(1, gene_symbol_col_name, gene_symbol)
            if len(gene_expression) < 1:
                gene_expression = df_donor_data_gene_probes.copy()
            else:
                gene_expression = gene_expression.merge(df_donor_data_gene_probes, how='outer')
                # gene_expression = gene_expression.assign(expression = gene_expression[df_gene[probe_col_name]].mean(axis=1))
    gene_expression = gene_expression.assign(expression=gene_expression[df_probe_ids[probe_col_name]].mean(axis=1))
    return gene_expression


def ai_extract(gene_symbols,df_donor_data, df_probes, df_wells2MNI):
    #["NTRK2", "MAG"]
    # given genes, donor data, probe mapping, and well mapping to MNI space, return dataframe of gene expression and coordinates
    genes2probes = get_probe_ids(gene_symbols,df_probes)
    gene_expression = get_gene_expression_multi(genes2probes,df_donor_data)
    gene_expression_coords = get_MNI_coords(gene_expression,df_wells2MNI)
    return gene_expression_coords
