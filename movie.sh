#!/bin/bash

ffmpeg -r 10 -i tl_psi_t/Rplot%03d.png -vcodec h264 -y -pix_fmt yuv420p psi_phi.mp4

ffmpeg -r 10 -i tl_psi_eta/Rplot%03d.png -vcodec h264 -y -pix_fmt yuv420p psi_eta.mp4

ffmpeg -y -i psi_eta.mp4 -i psi_phi.mp4 -filter_complex "[0:v]setpts=PTS-STARTPTS, pad=iw*2:ih[bg]; [1:v]setpts=PTS-STARTPTS[fg]; [bg][fg]overlay=w" psi_eta_phi.mp4