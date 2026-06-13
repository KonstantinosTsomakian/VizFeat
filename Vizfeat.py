import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import pyBigWig
import matplotlib.patches as patches





# Define a function that plots features from a given bed file

def plot_features_from_bed(bed_df, region, ax, parameter_list = [1,1 * (1 * 4), 1, 0, 0], name = ''):
    # bed_df is a pandas dataframe where the first 3 columns are chr, start, end
    # region is the exact same as the rest of the functions
    # ax is the matplolib ax where the features will e ploted
    # parameter list is the list of parameters for the dimensions of the plot
    scalefactor, step, width_scalefactor, y, total_step = parameter_list 
    
    width = 0.2 * width_scalefactor
    staring_point_on_y_axis = total_step + (width * 2)
    features_to_plot = bed_df[(bed_df.chr == region['chr'])]
    features_to_plot = features_to_plot[((features_to_plot.start >= region['start']) & (features_to_plot.start <= region['end'])) | ((features_to_plot.end >= region['start']) & (features_to_plot.end <= region['end']))]
    n_of_features = features_to_plot.shape[0]

    region_start = region['start']
    region_end = region['end']
    total_step+=10
    for i in range(n_of_features):
        
        start = features_to_plot.iloc[i,1] / 1000
        end = features_to_plot.iloc[i,2]  / 1000
        # print(start, end)
        
        # rect = patches.Rectangle((start, ((y + 1) * step) + 5), end - start, width, fill=True, edgecolor='black', linewidth=2, color = 'purple')
        rect = patches.Rectangle((start, total_step), end - start, width, fill=True, edgecolor='black', linewidth=2, color = 'purple')

        ax.add_patch(rect)
    finishing_point_on_y_axis = total_step + (width) + 5 # +5 because when plotting the hison mods i have an initial spacing of 10
    ax.vlines(x= (region_start / 1000) -5, ymin = staring_point_on_y_axis, ymax=finishing_point_on_y_axis, linewidth = 5, color = 'black')
    ax.text(s = name,x=(region_start / 1000) -10, y = staring_point_on_y_axis + ((finishing_point_on_y_axis - staring_point_on_y_axis) / 2), va = 'center', ha = 'right')
    parameter_list = [scalefactor, step, width_scalefactor, y, total_step ]
    print(f'finishing at {finishing_point_on_y_axis, staring_point_on_y_axis}')
    return(ax, parameter_list)



# Define a function to plot transcripts

def plot_transcripts(ax, region, transcript_annotation_df, scalefactor = 1,width_scalefactor=1,n_beds=0, total_step = 0, name = ''):
    staring_point_on_y_axis = total_step
    #transcript annotation df is a pandas dataframe where the columns are ['chr', 'type', 'start', 'end', 'strand', 'gene_id', 'transcript_id']
    # type refers to gene, exon, transcript
    # region is a dictionary with entries like {'chr': '15', 'start': 92759884, 'end': 92859884}
    # scalefactor adjusts the width of the transcript line
    # width scalefactor adjusts the height of the bars that define the exon regions
    # transcript_annotation_df is a pandas dataframe that stores chromosome name, the type, the start, the end, the strand, the gene id and the transcript id.
    # For each transcript there is an entry for a transcript and an entry for each one the exons that the transcript has.
    # Total_step is the height of the track
    


    transcript_annotation_df.chr = transcript_annotation_df.chr.astype(str)
    features_to_plot = transcript_annotation_df[(transcript_annotation_df.chr == str(region['chr']))]
    features_to_plot = features_to_plot[((features_to_plot.start >= region['start']) & (features_to_plot.start <= region['end'])) | ((features_to_plot.end >= region['start']) & (features_to_plot.end <= region['end']))]
    # print(features_to_plot)
    total_transcripts_to_plot = features_to_plot[features_to_plot.type == 'transcript'].transcript_id.to_list()
    # print(features_to_plot)
    n_of_transcripts = len(list(set(total_transcripts_to_plot)))
    print(f'{n_of_transcripts} transcripts found in the requested region.')

    # fig, ax = plt.subplots(figsize = (22, 11))

    region_start = region['start']
    region_end = region['end']
    

    step = 1 * (scalefactor * 4)
    # ax.set_ylim(ymin= -0.5, ymax  = (step * n_of_transcripts) + (n_beds+1))
    # ax.set_xlim(xmin = region_start  / 1000	,xmax = region_end  / 1000)
    chromosome = region['chr']
    

    
    for y, transcript_to_plot in enumerate(total_transcripts_to_plot):
        
        # print(transcript_to_plot)
        

        transcript_df = features_to_plot[features_to_plot.transcript_id == transcript_to_plot]
        transcript_start = transcript_df[transcript_df.type == 'transcript'].start.item() / 1000
        transcript_end = transcript_df[transcript_df.type == 'transcript'].end.item() / 1000
        transcript_strand = transcript_df[transcript_df.type == 'transcript'].strand.item()
        transcript_len = transcript_end - transcript_start
        # print(transcript_len)


        


        for i in range(transcript_df.shape[0]):
            width = 0.2 * width_scalefactor
            if transcript_df.iloc[i,1] == 'exon':
                start = transcript_df.iloc[i,2] / 1000
                end = transcript_df.iloc[i,3]  / 1000
                # print(start, end, y * step)
                rect = patches.Rectangle((start, y * step), end - start, width, fill=True, edgecolor='black', linewidth=2, color = 'black')
                ax.add_patch(rect)
                total_step = y * step

        line_xmin = region_start / 1000 if transcript_start < region_start / 1000 else transcript_start
        line_xmax = region_end / 1000 if transcript_end > region_end / 1000 else transcript_end
        ax.hlines(y = (y * step) + (width / 2), xmin=line_xmin, xmax = line_xmax, color = 'black', linewidth = 0.8 / scalefactor)


        if transcript_len >= 15:
            n_arrows = 10
        else:
            n_arrows = 5
        positions = np.linspace(line_xmin, line_xmax, n_arrows)  # 20 arrows
        # print(transcript_strand)
        if transcript_strand == '+':
            positions = positions[1:]
            arrowstyle = "->"
        elif transcript_strand == '-':
            positions = positions[:-1]
            arrowstyle = "<-"


        for x in positions[1:]:
            
            ax.annotate(
                "",
                xy=(x+0.00000000001, (y * step) + (width / 2)),   # arrow end
                xytext=(x, (y * step) + (width / 2)),            # arrow start
                arrowprops=dict(
                    arrowstyle=arrowstyle,
                    color="black",
                    linewidth=0.5
                    )
                )

        label_position = transcript_end + 1  # in kb if axis is in kb
        label_position = min(label_position, (region_end / 1000) + 1)  # leave small margin
        # print(label_position, region_end)
        ax.text(s = transcript_to_plot, x = label_position,y= (y * step) + width / 2,
                va='center', fontsize = 10)
    # ax.set_title(f'chr{chromosome}:{region_start}-{region_end}\n{(region_end - region_start) // 1000} kb',
    #         fontsize = 20)


    for spine in ax.spines.values():
        spine.set_visible(False)

    # ax.set_yticks([])
    
    
    # plt.close(fig)
    if n_of_transcripts == 0:
        y = 0
    finishing_point_on_y_axis = total_step + (width)
    ax.vlines(x= (region_start / 1000) -5, ymin = staring_point_on_y_axis, ymax=finishing_point_on_y_axis, linewidth = 5, color = 'black')
    ax.text(s = name,x=(region_start / 1000) -10, y = staring_point_on_y_axis + ((finishing_point_on_y_axis - staring_point_on_y_axis) / 2), va = 'center', ha = 'right')

    parameter_list=[scalefactor, step, width_scalefactor, y, total_step]
    # step is the step to move upwards in the y axis to plot each transcript
    # scalefactor controls the step and the width of the line that connects the exons. the bigger the scalefactor the bigger the step and the smaller theline
    # width scalefactor controls the width of the exon
    # y is the the y of the last plotted transcript and it is used to know where to plot any additional features
    return(ax, parameter_list)



# Define a function that takes as input a bw file and returns a numpy array of a the desired genomic region

def compute_binned_signal(bw,chrom,start,end,n_bins):
    # This function takes as input a bigwig file and returns the signal in bins as a numpy array.
    # chrom, start, end correspond to the chromosome to coordinates
    # n_bins is the desired number of bins

    import numpy
    # Compute bin edges
    bins = np.linspace(start, end, n_bins + 1, dtype=int)

    # Compute average signal in each bin
    binned_signal = []
    for i in range(n_bins):
        bin_start = bins[i]
        bin_end = bins[i+1]
        
        # pyBigWig values returns a list of values; missing values are None
        values = bw.values(chrom, bin_start, bin_end, numpy=True)
        
        # Convert None to 0 if desired
        values = np.nan_to_num(values)
        
        # Compute mean signal in this bin
        binned_signal.append(values.mean())

    binned_signal = np.array(binned_signal)
    return(binned_signal)


# Define a function that plots the genomic signal of the given genomic region

def add_hist_signal(ax, hist_signal_in_region,region, parameter_list = [1,1 * (1 * 4), 1, 0,0], color = 'red', name = ''):


	# This function takes as input the binned signal and returns the visualization
    # ax is the matplotlib ax
    # hist_signal_in_region is the the numpy array of the signal in the desired region
    # region is the dictionary that stores the chromosome the start and the end of the region
    # parameter list

    scalefactor, step, width_scalefactor, y, total_step = parameter_list
    staring_point_on_y_axis = total_step +9
    region_start = region['start']
    region_end = region['end']

    # print(max(hist_signal_in_region / 3))
    total_step+=10
    

    
    x = np.linspace(region_start / 1000, region_end / 1000, len(hist_signal_in_region))
    ax.fill_between(x = x, y2 = (hist_signal_in_region / 3) + total_step, y1 = total_step, color=color, alpha=0.5, label = name)
    total_step = max(total_step, total_step+max(hist_signal_in_region / 3).item())
    finishing_point_on_y_axis = total_step
    ax.vlines(x= (region_start / 1000) -5, ymin = staring_point_on_y_axis, ymax=finishing_point_on_y_axis + 1, linewidth = 5, color = 'black')
    ax.text(s = name,x=(region_start / 1000) -10, y = staring_point_on_y_axis + ((finishing_point_on_y_axis - staring_point_on_y_axis) / 2), va = 'center', ha = 'right')
    parameter_list = [scalefactor, step, width_scalefactor, y, total_step]
    return(ax, parameter_list)


