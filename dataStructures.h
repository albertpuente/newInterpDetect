#include <algorithm> // sort, max, min
#include <vector>
#include <climits>
#include <queue>
using namespace std;

struct Spike {
	int t;
	int chX, chY;
	unsigned short v;
    bool active;
	float distance(int chX2, int chY2) {
		return sqrt(pow((float) chX - (float) chX2, 2) + 
					pow((float) chY - (float) chX2, 2)   );
	}
};

struct Entry {
	int t;
	int x, y;
};

struct CheckSpace {
	queue<Entry> entries;
	vector<Spike>** spikes;
	int xSize, ySize;
	int maxD, maxT;
	int chunkSize;

	void initialise(float d, int t, int xS, int yS, int chunkS) {
		maxD = d;
		maxT = t;
		xSize = xS;
		ySize = yS;
		chunkSize = chunkS; // Must divide xSize and ySize

		spikes = new vector<Spike>*[xS/chunkSize];
		for (int i = 0; i < xS/chunkSize; i++) {
			spikes[i] = new vector<Spike>[yS/chunkSize];
			for (int j = 0; j < yS/chunkSize; j++) {
				spikes[i][j] = vector<Spike>();
			}
		}		
	}

	void addSpike(int t, int chX, int chY, bool status) {
		Spike s = {t, chX, chY, status};
		int x = chX/chunkSize;
		int y = chY/chunkSize;
		spikes[x][y].push_back(s);
		
		Entry e = {t, x, y};
		entries.push(e);
	}

	void pruneOldSpikes(int t) {
		while (entries.front().t < t - maxT) {
			Entry e = entries.front();
			auto it = spikes[e.x][e.y].begin();
			while (it != spikes[e.x][e.y].end()) {
				if ((*it).t < t - maxT)
                    // Print to file before deleting it
                    // TO-DO
					it = spikes[e.x][e.y].erase(it);
				else 
					it++;
			}			
			entries.pop();
		}
	}

	bool collides(int chX, int chY) {
		// This function checks if there has been a spike in the region. A maximum of 9 
		// chunks will be checked (less when the channel is next to a boundary).		

		// Chunk where the channel lies
		int x = chX/chunkSize;
		int y = chY/chunkSize;

		// In addition, check (if needed) the surrounding 8 chunks:
		// O O O
		// O X O 
		// O O O
        bool collision = false;
		for (int i = max(x - 1, 0); i < min(x + 1, xSize); i++) {
			for (int j = max(y - 1, 0); j < min(y + 1, ySize); j++) {				
				// Check all Spikes in chunk
				for (Spike & s : spikes[i][j]) {
					if (s.distance(chX, chY) < maxD) { 
                        collision = true;
						// Dactivate smaller spikes
                        // TO-DO
                    }
				}
			}
		}
		return false;
	}
};