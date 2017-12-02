#!/bin/bash 

for i in *.avi ; do 
	mp4_file=$(basename "${i/.avi}").mp4
	if [ -f "$mp4_file" ]
	then
		echo "$mp4_file created already"
	else	
		ffmpeg -i "$i" -c:v libx264 -crf 19 -preset slow -c:a libfdk_aac -b:a 192k -ac 2  "$mp4_file"
	fi
done


