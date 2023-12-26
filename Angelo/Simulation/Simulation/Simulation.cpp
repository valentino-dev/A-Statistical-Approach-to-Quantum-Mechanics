// Simulation.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>

using namespace std;

constexpr size_t TIME_COUNT = 1000;
constexpr size_t SITES_COUNT = 1000;

constexpr double MU = 1;
constexpr double LAMBDA = 1;

constexpr double A = 1;
constexpr double M_0 = 1;
constexpr double Delta = 1;
constexpr int n = 10;

double action(double positions[SITES_COUNT])
{
    double temp_sum = 0;
    for (size_t i = 0; i < SITES_COUNT-1; i++)
        temp_sum += M_0 / 2 * pow(positions[i] - positions[i + 1], 2) / A + (pow(MU, 2) / 2 * pow(positions[i], 2) + LAMBDA * pow(positions[i], 4)) * A;

    return temp_sum;
}

int main()
{
    double crystal[TIME_COUNT][SITES_COUNT]{ 0 };

    for (size_t j = 0; j < TIME_COUNT; j++) {
        for (size_t _ = 0; _ < n; _++) {
            double new_probabilities[SITES_COUNT]{ 0 }; // braucht noch eine uniform probablity distribution

            double d_action = action(new_probabilities) - action(crystal[j]);
            if (d_action < 0)
                for (size_t i = 0; i < SITES_COUNT; i++)
                    crystal[j][i] = new_probabilities[i];
            

            bool alles_im_Ramen = true;
            for (size_t i = 0; i < SITES_COUNT; i++)
                if (!(crystal[j][i] - Delta <= new_probabilities[i] <= crystal[j][i] + Delta)) alles_im_Ramen = false;

            if (alles_im_Ramen)
                for (size_t i = 0; i < SITES_COUNT; i++)
                    crystal[j][i] = new_probabilities[i];
        }
    }
}

// Programm ausführen: STRG+F5 oder Menüeintrag "Debuggen" > "Starten ohne Debuggen starten"
// Programm debuggen: F5 oder "Debuggen" > Menü "Debuggen starten"

// Tipps für den Einstieg: 
//   1. Verwenden Sie das Projektmappen-Explorer-Fenster zum Hinzufügen/Verwalten von Dateien.
//   2. Verwenden Sie das Team Explorer-Fenster zum Herstellen einer Verbindung mit der Quellcodeverwaltung.
//   3. Verwenden Sie das Ausgabefenster, um die Buildausgabe und andere Nachrichten anzuzeigen.
//   4. Verwenden Sie das Fenster "Fehlerliste", um Fehler anzuzeigen.
//   5. Wechseln Sie zu "Projekt" > "Neues Element hinzufügen", um neue Codedateien zu erstellen, bzw. zu "Projekt" > "Vorhandenes Element hinzufügen", um dem Projekt vorhandene Codedateien hinzuzufügen.
//   6. Um dieses Projekt später erneut zu öffnen, wechseln Sie zu "Datei" > "Öffnen" > "Projekt", und wählen Sie die SLN-Datei aus.
