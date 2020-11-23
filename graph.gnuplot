set terminal png font "Miso, 16"
set output 'graph.png'

set style data histogram
set style histogram cluster gap 1

set style fill solid border rgb "#ffffff"
set auto x
set title "Throughput (packet size=1280) Intel Core i5-8400 @ 2.8GHz"
set xlabel "K"
set ylabel "Mb/s"
set yrange [0:*]
plot 'graph.dat' using 2:xtic(1) lc rgb "#00ccff" title  col, \
        '' using 3:xtic(1) lc rgb "#00cccc" title col, \
        '' using 4:xtic(1) lc rgb "#0066ff" title col, \
        '' using 5:xtic(1) lc rgb "#004586" title col
