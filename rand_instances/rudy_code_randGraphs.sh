# randomly generated graphs

# density between [0 100]
# seed as 0

# density 20 %, n \in {100, 200, 300, 400,
# 500, 600, 700, 800, 900, 1000}

for i in {100..1000..100}
do
./rudy -rnd_graph "$i" 20 0 -random 0 100 0 >rndgraph_20_n$i.txt
done

# density 50 %, n \in {100, 200, 300, 400,
# 500, 600, 700, 800, 900, 1000}

for i in {100..1000..100}
do
./rudy -rnd_graph "$i" 50 0 -random 0 100 0 >rndgraph_50_n$i.txt
done

# density 80 %, n \in {100, 200, 300, 400,
# 500, 600, 700, 800, 900, 1000}

for i in {100..1000..100}
do
./rudy -rnd_graph "$i" 80 0 -random 0 100 0 >rndgraph_80_n$i.txt
done
