#include <SFML/Graphics.hpp>
#include <SFML/Graphics/Text.hpp>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <set>
#include <string>
#include <map>
#include <sstream>

void drawArrow(sf::RenderWindow &window, const sf::Font& font, const sf::Vector2f &start, const sf::Vector2f &end,
               const sf::Color &color, const std::string& weightStr) {
    sf::Vector2f dir = end - start;
    float fullLen = std::sqrt(dir.x*dir.x + dir.y*dir.y);
    if (fullLen <= 0.0001f) return;
    sf::Vector2f unit = dir / fullLen;

    float nodeOffset = 22.f;
    float arrowHeadSize = 18.f;
    float thickness = 5.f;

    sf::Vector2f newStart = start + unit * nodeOffset;
    sf::Vector2f newEnd = end - unit * nodeOffset;

    float bodyLen = std::sqrt((newEnd.x - newStart.x)*(newEnd.x - newStart.x) + (newEnd.y - newStart.y)*(newEnd.y - newStart.y));
    if (bodyLen <= 0.0001f) return;

    sf::RectangleShape body;
    body.setSize({ bodyLen, thickness });
    body.setOrigin({ 0.f, thickness/2.f });
    body.setPosition(newStart);

    float angleDeg = std::atan2(unit.y, unit.x) * 180.f / 3.14159265f;

    body.setRotation(sf::degrees(angleDeg));
    body.setFillColor(color);
    window.draw(body);

    sf::ConvexShape head;
    head.setPointCount(3);
    head.setPoint(0, sf::Vector2f(0.f, 0.f));
    head.setPoint(1, sf::Vector2f(-arrowHeadSize, arrowHeadSize/2.f));
    head.setPoint(2, sf::Vector2f(-arrowHeadSize, -arrowHeadSize/2.f));
    head.setFillColor(color);
    head.setPosition(newEnd);
    head.setRotation(sf::degrees(angleDeg));
    window.draw(head);

    sf::Text weightText(font);
    weightText.setString(weightStr);
    weightText.setCharacterSize(16);
    weightText.setFillColor(sf::Color::White);

    sf::Vector2f textPos = newStart + unit * (bodyLen / 2.f);
    sf::Vector2f perpendicular = sf::Vector2f(-unit.y, unit.x);
    textPos += perpendicular * 15.f;

    sf::FloatRect textRect = weightText.getLocalBounds();
    weightText.setOrigin({textRect.position.x + textRect.size.x / 2.0f,
                         textRect.position.y + textRect.size.y / 2.0f});
    weightText.setPosition(textPos);
    window.draw(weightText);
}

std::vector<int> bfsOrder(const std::vector<std::vector<int>>& adj, std::vector<int>& parent, int start = 0) {
    int n = (int)adj.size();
    std::vector<int> vis(n, 0), ord;
    std::queue<int> q;
    q.push(start);
    vis[start] = 1;
    parent[start] = -1;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        ord.push_back(u);
        for (int v = 0; v < n; ++v) {
            if (adj[u][v] && !vis[v]) {
                vis[v] = 1;
                parent[v] = u;
                q.push(v);
            }
        }
    }
    return ord;
}

std::vector<int> dijkstraOrder(const std::vector<std::vector<int>>& adj, std::vector<int>& parent,
                               std::vector<long long>& dist, int start = 0) {
    int n = (int)adj.size();
    dist.assign(n, LLONG_MAX);
    std::vector<int> vis(n, 0), ord;
    parent.assign(n, -1);

    dist[start] = 0;

    for (int i = 0; i < n; ++i) {
        int u = -1;
        for (int j = 0; j < n; ++j) {
            if (!vis[j] && (u == -1 || dist[j] < dist[u])) {
                u = j;
            }
        }

        if (u == -1 || dist[u] == LLONG_MAX) break;

        vis[u] = 1;
        ord.push_back(u);

        for (int v = 0; v < n; ++v) {
            if (adj[u][v] != 0 && dist[u] + adj[u][v] < dist[v]) {
                dist[v] = dist[u] + adj[u][v];
                parent[v] = u;
            }
        }
    }
    return ord;
}

bool bellmanFord(const std::vector<std::vector<int>>& adj, std::vector<int>& parent,
                   std::vector<long long>& dist, int start = 0) {
    int n = (int)adj.size();
    dist.assign(n, LLONG_MAX);
    parent.assign(n, -1);

    dist[start] = 0;

    for (int i = 1; i <= n - 1; ++i) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                int weight = adj[u][v];
                if (weight != 0 && dist[u] != LLONG_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                }
            }
        }
    }

    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            int weight = adj[u][v];
            if (weight != 0 && dist[u] != LLONG_MAX && dist[u] + weight < dist[v]) {
                return true;
            }
        }
    }

    return false;
}

int minKey(int n, const std::vector<int>& key, const std::vector<bool>& inMST) {
    int minVal = INT_MAX, minIndex = -1;
    for (int v = 0; v < n; ++v) {
        if (!inMST[v] && key[v] < minVal) {
            minVal = key[v];
            minIndex = v;
        }
    }
    return minIndex;
}

std::vector<int> mstPrimsOrder(const std::vector<std::vector<int>>& adj, std::vector<int>& parent, int start = 0) {
    int n = (int)adj.size();
    std::vector<int> key(n, INT_MAX);
    std::vector<bool> inMST(n, false);
    std::vector<int> ord;

    key[start] = 0;
    parent[start] = -1;

    for (int count = 0; count < n; ++count) {
        int u = minKey(n, key, inMST);
        if (u == -1) break;

        inMST[u] = true;
        ord.push_back(u);

        for (int v = 0; v < n; ++v) {
            if (adj[u][v] > 0 && !inMST[v] && adj[u][v] < key[v]) {
                parent[v] = u;
                key[v] = adj[u][v];
            }
        }
    }
    return ord;
}

float heuristic(int nodeA, int nodeB, const std::vector<sf::Vector2f>& pos) {
    sf::Vector2f dir = pos[nodeB] - pos[nodeA];
    return std::sqrt(dir.x * dir.x + dir.y * dir.y);
}

std::vector<int> aStarOrder(const std::vector<std::vector<int>>& adj, std::vector<int>& parent,
                            std::vector<long long>& gScore, const std::vector<sf::Vector2f>& pos,
                            int start, int dest) {
    int n = (int)adj.size();
    gScore.assign(n, LLONG_MAX);
    std::vector<float> fScore(n, std::numeric_limits<float>::infinity());
    std::vector<int> vis(n, 0), ord;
    parent.assign(n, -1);

    gScore[start] = 0;
    fScore[start] = heuristic(start, dest, pos);

    for (int i = 0; i < n; ++i) {
        int u = -1;
        float minF = std::numeric_limits<float>::infinity();

        for (int j = 0; j < n; ++j) {
            if (!vis[j] && fScore[j] < minF) {
                minF = fScore[j];
                u = j;
            }
        }

        if (u == -1) break;

        vis[u] = 1;
        ord.push_back(u);

        if (u == dest) break;

        for (int v = 0; v < n; ++v) {
            if (adj[u][v] != 0) {
                long long tentative_gScore = gScore[u] + adj[u][v];
                if (tentative_gScore < gScore[v]) {
                    parent[v] = u;
                    gScore[v] = tentative_gScore;
                    fScore[v] = (float)gScore[v] + heuristic(v, dest, pos);
                }
            }
        }
    }
    return ord;
}

void visualizeGraph(const std::vector<std::vector<int>>& adj,
                    const std::vector<int>& order,
                    const std::vector<int>& parent,
                    const std::set<std::pair<int, int>>& shortestPathEdges,
                    const std::vector<long long>& dist,
                    const std::vector<sf::Vector2f>& pos,
                    const std::string& title,
                    int destIndex,
                    long long totalCost)
{
    int n = (int)adj.size();
    sf::RenderWindow window(sf::VideoMode(sf::Vector2u{800u, 600u}), title);

    sf::Font font;
    if (!font.openFromFile("C:/Windows/Fonts/arial.ttf")) {
        std::cerr << "Error: Could not load font 'arial.ttf'.\n";
        return;
    }

    float radius = 20.f;

    std::vector<sf::Color> nodeColor(n, sf::Color(200,200,200));
    sf::Color visitedNodeColor = sf::Color::Cyan;
    sf::Color visitingArrowColor = sf::Color::Yellow;
    sf::Color visitedArrowColor = sf::Color::Red;
    sf::Color finalPathColor = sf::Color::Green;

    sf::Text titleText(font);
    titleText.setString(title);
    titleText.setCharacterSize(24);
    titleText.setFillColor(sf::Color::White);
    titleText.setStyle(sf::Text::Bold);
    sf::FloatRect titleRect = titleText.getLocalBounds();
    titleText.setOrigin({titleRect.position.x + titleRect.size.x / 2.0f,
                         titleRect.position.y + titleRect.size.y / 2.0f});
    titleText.setPosition({400.f, 30.f});

    sf::Text infoText(font);
    infoText.setCharacterSize(20);
    infoText.setFillColor(sf::Color::White);
    infoText.setPosition({30.f, 70.f});

    sf::Text statusText(font);
    statusText.setCharacterSize(20);
    statusText.setFillColor(sf::Color::White);
    statusText.setPosition({30.f, 550.f});

    sf::Clock clk;
    size_t step = 0;
    bool animationComplete = false;

    while (window.isOpen()) {
        for (auto ev = window.pollEvent(); ev.has_value(); ev = window.pollEvent()) {
            if (ev->is<sf::Event::Closed>()) window.close();
        }

        if (step < order.size() && clk.getElapsedTime().asSeconds() > 0.5f) {
            nodeColor[order[step]] = visitedNodeColor;
            step++;
            clk.restart();
        }
        if (step == order.size()) {
            animationComplete = true;
        }

        if (animationComplete) {
            statusText.setString("Done! Final path in green.");
        } else {
            std::string status = "Step: " + std::to_string(step) + " / " + std::to_string(order.size());
            if (step < order.size()) {
                status += " (Visiting Node " + std::to_string(order[step] + 1) + ")";
            }
            statusText.setString(status);
        }

        if (title.find("MST") != std::string::npos) {
            infoText.setString("Total MST Cost: " + std::to_string(totalCost));
        } else if (title.find("BFS") != std::string::npos) {
            infoText.setString("Total Cost: N/A (unweighted)");
        } else if (title.find("Dijkstra") != std::string::npos || title.find("Bellman-Ford") != std::string::npos || title.find("A*") != std::string::npos) {
            if (destIndex != -1) {
                std::string dStr = (dist[destIndex] == LLONG_MAX) ? "inf" : std::to_string(dist[destIndex]);
                infoText.setString("Shortest Path to Node " + std::to_string(destIndex + 1) + ": " + dStr);
            } else {
                infoText.setString("Destination: Not Selected");
            }
        }

        window.clear(sf::Color(30,30,40));

        if (animationComplete) {
            for (int j = 0; j < n; ++j) {
                int i = parent[j];
                if (i != -1 && nodeColor[j] == visitedNodeColor) {
                    if (shortestPathEdges.find({i, j}) == shortestPathEdges.end()) {
                        drawArrow(window, font, pos[i], pos[j], visitedArrowColor, std::to_string(adj[i][j]));
                    }
                }
            }

            for (const auto& edge : shortestPathEdges) {
                int i = edge.first;
                int j = edge.second;
                drawArrow(window, font, pos[i], pos[j], finalPathColor, std::to_string(adj[i][j]));
            }

        } else {
            for (size_t k = 0; k < step; ++k) {
                int j = order[k];
                int i = parent[j];

                if (i != -1) {
                    sf::Color arrowColor;
                    if (k == step - 1) {
                        arrowColor = visitingArrowColor;
                    } else {
                        arrowColor = visitedArrowColor;
                    }
                    drawArrow(window, font, pos[i], pos[j], arrowColor, std::to_string(adj[i][j]));
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            sf::CircleShape node(radius);
            node.setOrigin({radius, radius});
            node.setPosition(pos[i]);
            node.setFillColor(nodeColor[i]);
            node.setOutlineThickness(2.f);
            node.setOutlineColor(sf::Color(240, 240, 240));
            window.draw(node);

            sf::Text text(font);
            text.setString(std::to_string(i + 1));
            text.setCharacterSize(18);
            text.setFillColor(sf::Color::Black);

            sf::FloatRect textRect = text.getLocalBounds();
            text.setOrigin({textRect.position.x + textRect.size.x / 2.0f,
                           textRect.position.y + textRect.size.y / 2.0f});
            text.setPosition(pos[i]);
            window.draw(text);

            if (title.find("Dijkstra") != std::string::npos || title.find("Bellman-Ford") != std::string::npos || title.find("A*") != std::string::npos) {
                sf::Text distText(font);
                distText.setCharacterSize(16);
                distText.setFillColor(sf::Color::White);

                std::string dStr = "d = ";
                if (dist[i] == LLONG_MAX) {
                    dStr += "inf";
                } else {
                    dStr += std::to_string(dist[i]);
                }

                distText.setString(dStr);

                sf::FloatRect distRect = distText.getLocalBounds();
                distText.setOrigin({distRect.position.x + distRect.size.x / 2.0f,
                                    distRect.position.y + distRect.size.y / 2.0f});
                distText.setPosition(pos[i] + sf::Vector2f(0.f, radius + 12.f));
                window.draw(distText);
            }
        }

        window.draw(titleText);
        window.draw(infoText);
        window.draw(statusText);
        window.display();
    }
}

void runGraphAlgorithms() {
    int n;
    std::cout << "Enter number of nodes (max 10 recommended): ";
    if (!(std::cin >> n) || n <= 0) {
        std::cerr << "Invalid number\n";
        return;
    }

    std::vector<std::vector<int>> graph(n, std::vector<int>(n));
    std::cout << "Enter adjacency matrix (" << n << "x" << n << "), 0 for no edge:\n";
    std::cout << "(Use negative numbers for Bellman-Ford)\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            std::cin >> graph[i][j];

    std::vector<sf::Vector2f> pos(n);
    float cx = 400.f, cy = 300.f, R = 200.f;
    for (int i = 0; i < n; ++i) {
        float ang = i * 2.f * 3.14159265f / std::max(1, n);
        pos[i] = sf::Vector2f(cx + R * std::cos(ang), cy + R * std::sin(ang));
    }

    std::cout << "Choose algorithm:\n1. Dijkstra\n2. BFS\n3. Bellman-Ford\n4. MST (Prim's)\n5. A* Search\n> ";
    int ch; std::cin >> ch;

    std::vector<int> parent(n, -1);
    std::vector<int> order;
    std::string title;
    bool isShortestPathAlgo = false;
    std::vector<long long> dist(n, LLONG_MAX);
    long long totalCost = 0;
    int destIndex = -1;

    if (ch == 1) {
        order = dijkstraOrder(graph, parent, dist, 0);
        title = "Dijkstra Visualization";
        isShortestPathAlgo = true;
    } else if (ch == 2) {
        order = bfsOrder(graph, parent, 0);
        title = "BFS Visualization";
        isShortestPathAlgo = true;
    } else if (ch == 3) {
        bool hasNegativeCycle = bellmanFord(graph, parent, dist, 0);

        if (hasNegativeCycle) {
            std::cout << "\n\n*** WARNING: Negative weight cycle detected! ***\n";
            std::cout << "Shortest paths are not well-defined.\n\n";
        }

        for (int i = 0; i < n; ++i) {
            if (dist[i] != LLONG_MAX) {
                order.push_back(i);
            }
        }
        title = "Bellman-Ford Visualization";
        isShortestPathAlgo = true;
    } else if (ch == 4) {
        order = mstPrimsOrder(graph, parent, 0);
        title = "Prim's MST Visualization";
        isShortestPathAlgo = false;
        for (int i = 0; i < n; ++i) {
            if (parent[i] != -1) {
                totalCost += graph[parent[i]][i];
            }
        }
    } else if (ch == 5) {
        title = "A* Search Visualization";
        isShortestPathAlgo = true;
        std::cout << "A* Search requires a destination.\n";
        std::cout << "Enter destination node (1-" << n << "): ";
        int destNode;
        if (!(std::cin >> destNode) || destNode < 1 || destNode > n) {
            std::cerr << "Invalid destination node.\n";
            return;
        }
        destIndex = destNode - 1;
        order = aStarOrder(graph, parent, dist, pos, 0, destIndex);
    } else {
        std::cerr << "Invalid choice\n";
        return;
    }


    std::set<std::pair<int, int>> shortestPathEdges;

    if (isShortestPathAlgo && destIndex == -1) {
        std::cout << "Enter destination node (1-" << n << ") to highlight path (0 for none): ";
        int destNode;
        std::cin >> destNode;
        destIndex = destNode - 1;
    }

    if (isShortestPathAlgo && destIndex >= 0) {
        if (destIndex < n && (parent[destIndex] != -1 || destIndex == 0) ) {
            int curr = destIndex;
            while (parent[curr] != -1) {
                shortestPathEdges.insert({parent[curr], curr});
                curr = parent[curr];
            }
        } else if (destIndex != 0) {
            std::cout << "Node " << (destIndex + 1) << " is not reachable from node 1.\n";
        }
    }


    std::cout << "\nTraversal/Visit order: ";
    for (size_t i = 0; i < order.size(); ++i) {
        std::cout << (order[i] + 1);
        if (i + 1 < order.size()) std::cout << " -> ";
    }
    std::cout << "\n\nVisited nodes (step by step):\n";
    for (auto v : order) std::cout << "Node " << (v+1) << "\n";

    visualizeGraph(graph, order, parent, shortestPathEdges, dist, pos, title, destIndex, totalCost);
}

struct HuffmanNode {
    char data;
    unsigned freq;
    HuffmanNode *left, *right;

    HuffmanNode(char data, unsigned freq) : data(data), freq(freq), left(nullptr), right(nullptr) {}
};

struct compare {
    bool operator()(HuffmanNode* l, HuffmanNode* r) {
        return (l->freq > r->freq);
    }
};

HuffmanNode* buildHuffmanTree(std::map<char, unsigned>& freqMap) {
    std::priority_queue<HuffmanNode*, std::vector<HuffmanNode*>, compare> minHeap;

    for (auto pair : freqMap) {
        minHeap.push(new HuffmanNode(pair.first, pair.second));
    }

    while (minHeap.size() > 1) {
        HuffmanNode* left = minHeap.top(); minHeap.pop();
        HuffmanNode* right = minHeap.top(); minHeap.pop();

        HuffmanNode* top = new HuffmanNode('$', left->freq + right->freq);
        top->left = left;
        top->right = right;
        minHeap.push(top);
    }
    return minHeap.top();
}

void generateCodes(HuffmanNode* root, std::string str, std::map<char, std::string>& codes) {
    if (!root) return;

    if (root->data != '$') {
        codes[root->data] = str;
    }
    generateCodes(root->left, str + "0", codes);
    generateCodes(root->right, str + "1", codes);
}

std::string decodeHuffman(HuffmanNode* root, std::string s) {
    std::string ans = "";
    HuffmanNode* curr = root;
    for (int i=0; i<s.length(); i++) {
        if (s[i] == '0')
           curr = curr->left;
        else
           curr = curr->right;

        if (curr->left == nullptr && curr->right == nullptr) {
            ans += curr->data;
            curr = root;
        }
    }
    return ans;
}

void deleteTree(HuffmanNode* node) {
    if (node == nullptr) return;
    deleteTree(node->left);
    deleteTree(node->right);
    delete node;
}

void drawLine(sf::RenderWindow& window, sf::Vector2f p1, sf::Vector2f p2) {
    sf::Vertex line[] = { sf::Vertex{p1, sf::Color::White}, sf::Vertex{p2, sf::Color::White} };
    window.draw(line, 2, sf::PrimitiveType::Lines);
}

void drawTreeNodes(sf::RenderWindow& window, const sf::Font& font, HuffmanNode* node, float x, float y, float hSpacing) {
    if (node == nullptr) return;

    sf::CircleShape circle(20.f);
    circle.setOrigin({20.f, 20.f});
    circle.setPosition({x, y});
    circle.setFillColor(sf::Color::Cyan);
    circle.setOutlineColor(sf::Color::White);
    circle.setOutlineThickness(2.f);
    window.draw(circle);

    sf::Text text(font);
    std::string label;
    if (node->data == '$') {
        label = std::to_string(node->freq);
    } else if (node->data == ' ') {
        label = "' '";
    } else if (node->data == '\n') {
        label = "\\n";
    } else {
        label = std::string(1, node->data);
    }
    text.setString(label);
    text.setCharacterSize(18);
    text.setFillColor(sf::Color::Black);
    sf::FloatRect textRect = text.getLocalBounds();
    text.setOrigin({textRect.position.x + textRect.size.x / 2.0f,
                   textRect.position.y + textRect.size.y / 2.0f});
    text.setPosition({x, y});
    window.draw(text);

    float vSpacing = 80.f;
    if (node->left) {
        sf::Vector2f leftPos(x - hSpacing, y + vSpacing);
        drawLine(window, {x, y + 20.f}, {leftPos.x, leftPos.y - 20.f});
        drawTreeNodes(window, font, node->left, leftPos.x, leftPos.y, hSpacing / 2.f);
    }
    if (node->right) {
        sf::Vector2f rightPos(x + hSpacing, y + vSpacing);
        drawLine(window, {x, y + 20.f}, {rightPos.x, rightPos.y - 20.f});
        drawTreeNodes(window, font, node->right, rightPos.x, rightPos.y, hSpacing / 2.f);
    }
}

void drawCodeTable(sf::RenderWindow& window, const sf::Font& font, std::map<char, std::string>& codes) {
    std::stringstream ss;
    ss << "Huffman Codes:\n";
    for (auto pair : codes) {
        if (pair.first == ' ') ss << "' '";
        else if (pair.first == '\n') ss << "\\n";
        else ss << " " << pair.first;
        ss << " : " << pair.second << "\n";
    }

    sf::Text tableText(font);
    tableText.setString(ss.str());
    tableText.setCharacterSize(18);
    tableText.setFillColor(sf::Color::White);
    tableText.setPosition({30.f, 30.f});
    window.draw(tableText);
}

void visualizeHuffman(HuffmanNode* root, std::map<char, std::string>& codes) {
    sf::RenderWindow window(sf::VideoMode(sf::Vector2u{1200u, 800u}), "Huffman Tree Visualization");

    sf::Font font;
    if (!font.openFromFile("C:/Windows/Fonts/arial.ttf")) {
        std::cerr << "Error: Could not load font 'arial.ttf'.\n";
        return;
    }

    while (window.isOpen()) {
        for (auto ev = window.pollEvent(); ev.has_value(); ev = window.pollEvent()) {
            if (ev->is<sf::Event::Closed>()) window.close();
        }

        window.clear(sf::Color(30,30,40));
        drawTreeNodes(window, font, root, 600.f, 100.f, 300.f);
        drawCodeTable(window, font, codes);
        window.display();
    }
}

void runHuffman() {
    std::cout << "Enter a string to encode: ";
    std::string text;
    std::getline(std::cin >> std::ws, text);

    std::map<char, unsigned> freqMap;
    for (char c : text) {
        freqMap[c]++;
    }

    HuffmanNode* root = buildHuffmanTree(freqMap);

    std::map<char, std::string> codes;
    generateCodes(root, "", codes);

    std::cout << "\nHuffman Codes:\n";
    std::string encodedString = "";
    for (char c : text) {
        std::string codeStr = codes[c];
        encodedString += codeStr;

        if (c == ' ') std::cout << "' '";
        else if (c == '\n') std::cout << "\\n";
        else std::cout << " " << c;
        std::cout << " : " << codeStr << "\n";
    }

    std::cout << "\nOriginal string (" << text.length() << " chars): " << text << "\n";
    std::cout << "Encoded string (" << encodedString.length() << " bits): " << encodedString << "\n";

    std::string decodedString = decodeHuffman(root, encodedString);
    std::cout << "Decoded string: " << decodedString << "\n";

    visualizeHuffman(root, codes);

    deleteTree(root);
}


int main() {
    std::cout << "Choose an application:\n";
    std::cout << "1. Graph Algorithm Visualizer\n";
    std::cout << "2. Huffman Encoding Visualizer\n";
    std::cout << "> ";
    int choice;
    std::cin >> choice;

    if (choice == 1) {
        runGraphAlgorithms();
    } else if (choice == 2) {
        runHuffman();
    } else {
        std::cerr << "Invalid choice.\n";
        return 1;
    }

    return 0;
} 