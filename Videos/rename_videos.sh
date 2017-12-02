#!/bin/bash 

# S1 video 
mp4_file="example_microendoscopic_data.mp4"
new_name="S1 Video.mp4"
if [ -f "$mp4_file" ]; then
	mv "$mp4_file" "$new_name"
else
	echo "$mp4_file does not exist"
fi

# S2 video 
mp4_file="background_comparison.mp4"
new_name="S2 Video.mp4"
if [ -f "$mp4_file" ]; then
	mv "$mp4_file" "$new_name"
else
	echo "$mp4_file does not exist"
fi

# S3 video 
mp4_file="sim_initializatoin.mp4"
new_name="S3 Video.mp4"
if [ -f "$mp4_file" ]; then
	mv "$mp4_file" "$new_name"
else
	echo "$mp4_file does not exist"
fi


# S4 video 
mp4_file="sim_snr1_decompose.mp4"
new_name="S4 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S5 video 
mp4_file="sim_snr6_decompose.mp4"
new_name="S5 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S6 video 
mp4_file="striatum_decompose.mp4"
new_name="S6 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S7 video 
mp4_file="PFC_decompose.mp4"
new_name="S7 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S8 video 
mp4_file="PFC_overlapping.mp4"
new_name="S8 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S9 video 
mp4_file="HIPPOCAMPUS_decompose.mp4"
new_name="S9 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

# S10 video 
mp4_file="intervention.mp4"
new_name="S10 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi


# S11 video 
mp4_file="BNST_decompose.mp4"
new_name="S11 Video.mp4"
if [ -f "$mp4_file" ]; then
        mv "$mp4_file" "$new_name"
else
        echo "$mp4_file does not exist"
fi

