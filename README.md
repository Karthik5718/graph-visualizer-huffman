<mark>Use C++20 and SFML 3.0.1</mark>

A single-file C++ + SFML project that visualizes multiple graph algorithms and Huffman encoding/decoding.
Includes real-time animations, node/edge highlighting, path reconstruction, and interactive controls.

All algorithms are animated step-by-step:

1. Dijkstra’s Algorithm

Shows distance relaxation

Highlights current node and shortest-path edges

Displays final shortest path tree

2. BFS

Level-order exploration

Queue visualization

Highlights frontier layers

3. Bellman–Ford

Shows relaxation for every iteration

Detects negative cycles

Full step-by-step transitions

4. MST (Prim’s Algorithm)

Highlights minimum edge selected at each step

Builds MST tree live

Useful for teaching greedy algorithms

5. A Search*

Grid-based or graph-based (your choice)

g(n), h(n), and f(n) visualization

Highlights open and closed lists

Shows final optimal path

🗜️ Huffman Encoding / Decoding

1. Includes a fully visualized Huffman workflow:

2. Generates frequency table

3. Builds min-heap visually

4. Shows node merges and final tree

5. Displays binary codes for each char

6. Encodes and decodes text

Animations for each step
Visual Features (SFML)

1. Node animations, edge highlighting

2. Smooth transitions

3. Color themes for visited, processing, and result states

4. Real-time rendering using SFML shapes, text, lines
