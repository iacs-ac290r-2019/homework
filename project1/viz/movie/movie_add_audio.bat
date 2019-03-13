ffmpeg -safe 0 -f concat -i movie_files.txt -c copy movie_video.mp4
ffmpeg -i movie_video.mp4 -i chopin_25_01.mp3 -c:v copy -shortest movie.mp4