#include <cassert>
#include <iostream>
#include <random>

// Game settings

const int SimulationCount = 100000;   // How many simulations to run

const int GameCount = 16000; // How many games to run in total
const int GameChunk = 1000;  // How many games between testing stopping condition.

// Patch settings (CRITICAL SECTION!)

const double PatchEloMean = -2.0;              // Average ELO result when testing a single patch
const double PatchEloStdDev = 3.0;             // Standard deviation of ELO result when testing a single patch
const double PatchEloAcceptanceCriteria = 2.5; // Result patch needs to reach to be "accepted"

// One liners and prototyps

double elo_to_winperc(double eloDelta) {  return 1.0 / (1.0 + pow(10.0, - eloDelta / 400.0)); }
double winperc_to_elo(double winPerc)  {  return -400 * log10((1 - winPerc) / winPerc); }
double run_one_simulation(double winPerc, int& stoppedAfterGames);

// Global variables

std::random_device rd;
std::mt19937 gen(rd());

// Early stopping condition (CRITICAL SECTION!)

bool early_stop(int wonGames, int playedGames)
{
   if (playedGames == 8000)
       return wonGames < 4000;
   else
       return false;
}

// Main function

int main()
{
    int passedWithoutStoppingCond = 0;
    int passedWithStoppingCond = 0;
    int gamesSaved = 0;

    assert(GameCount % GameChunk == 0);

    std::normal_distribution<> d(PatchEloMean, PatchEloStdDev);
 
    for (int i = 0; i < SimulationCount ; i++)
    {
        int stoppedAfterGames;
        double eloGain = winperc_to_elo(run_one_simulation(elo_to_winperc(d(gen)), stoppedAfterGames));

        if (eloGain > PatchEloAcceptanceCriteria)
            passedWithoutStoppingCond++;
        
        if (eloGain > PatchEloAcceptanceCriteria && !stoppedAfterGames)
            passedWithStoppingCond++;

        if (stoppedAfterGames)
            gamesSaved += GameCount - stoppedAfterGames;
    }

    std::cout << "Games saved : " << gamesSaved << " / " << SimulationCount * GameCount << " ( " << double(gamesSaved) / (SimulationCount * GameCount) * 100 << "% )\n";
    std::cout << "Patches passed (without stopping condition) : " << passedWithoutStoppingCond << " / " << SimulationCount << " ( " << double(passedWithoutStoppingCond) / SimulationCount * 100 << "% )\n";
    std::cout << "Patches passed (with stopping condition )   : " << passedWithStoppingCond << " / " << SimulationCount << " ( " << double(passedWithStoppingCond) / SimulationCount * 100 << "% )\n";
    std::cout << "Good Patch discard rate : " << passedWithoutStoppingCond - passedWithStoppingCond << " / " << passedWithoutStoppingCond << " ( " << double(passedWithoutStoppingCond - passedWithStoppingCond) / passedWithoutStoppingCond * 100 << "% )\n"; 

    return 0;
}

// IN:      winPerc           -- Patch Real Winning percentage for a single game
// OUT:     stoppedAfterGames -- If zero, simulation was not stopped. If non-zero, is set to how many games were played before simulation was stopped.
// RETURNS: winning percentage of the patch for the whole simulation.

double run_one_simulation(double winPerc, int& stoppedAfterGames)
{
    int iterationCount = GameCount / GameChunk;
    int playedGames = 0, wonGames = 0;

    stoppedAfterGames = 0;

    // Approximate binomial distribution with normal distribution: B(x; n, p) ~= N(np, np(1-p))
    std::normal_distribution<> d(GameChunk * winPerc, sqrt(GameChunk * winPerc * (1.0 - winPerc)));

    for (int i = 0; i < iterationCount ; i++)
    {
        playedGames += GameChunk;
        wonGames += std::round(d(gen));

        // Stop test early?
        if (!stoppedAfterGames && early_stop(wonGames, playedGames))
            stoppedAfterGames = playedGames;
    }

    return double(wonGames) / playedGames;
}
