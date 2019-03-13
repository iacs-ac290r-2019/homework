rem ffmpeg -r 12 -f image2 -s 1500x600 -i frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p temp.mp4
ffmpeg -r 12 -f image2 -s 1500x600 -i temp_frames_mesh/frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p movie/temperature.mp4
ffmpeg -r 12 -f image2 -s 1500x600 -i p_frames/frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p movie/pressure.mp4
ffmpeg -r 12 -f image2 -s 1500x600 -i ux_frames/frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p movie/ux.mp4
ffmpeg -r 12 -f image2 -s 1500x600 -i uy_frames/frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p movie/uy.mp4
ffmpeg -r 12 -f image2 -s 1500x600 -i streamline_frames/reg_frame%04d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p movie/streamline.mp4