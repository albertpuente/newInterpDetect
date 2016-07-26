#include <algorithm> // sort, max, min, fill
#include <vector>
#include <climits>
#include <queue>
#include <iostream>
#include <mutex>
using namespace std;

struct Spike {
	int t;
	int chX, chY;
	int amp;
    bool relevant;
	inline float distance(int chX2, int chY2) {
		if (chX == chX2 and chY == chY2)
			return 0;
		else 
			return sqrt(pow((float) chX - (float) chX2, 2) + 
						pow((float) chY - (float) chY2, 2)   );
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

	std::mutex global_enqueue;

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

	void addSpike(int amp, int t, int chX, int chY, bool status) {
		Spike s = {t, chX, chY, amp, status};
		int x = chX/chunkSize;
		int y = chY/chunkSize;
		
		spikes[x][y].push_back(s);
		
		Entry e = {t, x, y};

		// This is the only possible concurrent access
		global_enqueue.lock();
		entries.push(e);
		global_enqueue.unlock();
	}

	vector<Spike> pruneOldSpikes(int t, int* fourChInterp, int* fiveChInterp) {
		auto spikesToPrint = vector<Spike>();
		while (not entries.empty() && entries.front().t < t - maxT) {			
			Entry e = entries.front();
			auto it = spikes[e.x][e.y].begin();
			while (it != spikes[e.x][e.y].end()) {
				if ((*it).t < t - maxT) {
                    // Print to file before deleting it
                    // TO-DO
					Spike spk = *it;
					if (spk.relevant) {
						spikesToPrint.push_back(spk);
					}
					it = spikes[e.x][e.y].erase(it);
				}
				else {
					it++;
				}
			}			
			entries.pop();
		}
		return spikesToPrint;
	}

	bool collides(int amp, int chX, int chY) {
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
		for (int i = max(x - 1, 0); i < min(x + 2, xSize/chunkSize); i++) {
			for (int j = max(y - 1, 0); j < min(y + 2, ySize/chunkSize); j++) {				
				// Check all Spikes in chunk
				for (Spike & s : spikes[i][j]) {					
					if (s.distance(chX, chY) < maxD) {
						if (s.amp < amp) {
							s.relevant = false;		
						}
						else {
							collision = true;
						}                    					
                    }
				}
			}
		}
		return collision;
	}

	void purge() {
		while (not entries.empty()) {			
			Entry e = entries.front();
			auto it = spikes[e.x][e.y].begin();
			while (it != spikes[e.x][e.y].end()) {
				it = spikes[e.x][e.y].erase(it);
			}	
			entries.pop();
		}
	}
};