#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "TSystem.h"

using namespace std;

void DownloadRunData(string files_to_download_str = "") {
    
    // 1. Configurazione dei percorsi e ID di Drive
    string drive_id = "1Wm7TsoTfgXSlrE2OFXad0jpF6V2T1P6a";
    string target_folder = "/home/emanuele/Universita/Paprica 2.9/Dataset/Run_Data";
    
    // Creiamo la cartella di destinazione se non esiste già
    if (gSystem->AccessPathName(target_folder.c_str())) {
        cout << "La cartella di destinazione non esiste. Creazione in corso..." << endl;
        gSystem->mkdir(target_folder.c_str(), true);
    }
    
    // 2. Controllo degli input da riga di comando
    if (files_to_download_str.empty()) {
        cerr << "Errore: Non hai specificato nessun file da scaricare!" << endl;
        cerr << "Uso: root -l 'DownloadRunData.C(\"file1.txt,file2.txt\")'" << endl;
        return;
    }
    
    // Parsing dei file separati da virgola
    vector<string> files_list;
    stringstream ss(files_to_download_str);
    string token;
    while (getline(ss, token, ',')) {
        // Rimuove eventuali spazi bianchi iniziali o finali accidentalmente inseriti
        token.erase(0, token.find_first_not_of(" "));
        token.erase(token.find_last_not_of(" ") + 1);
        if (!token.empty()) {
            files_list.push_back(token);
        }
    }
    
    cout << "\n=== Inizio Download Selettivo da Google Drive ===" << endl;
    cout << "Destinazione: " << target_folder << "\n" << endl;
    
    // 3. Loop di download tramite rclone
    int download_counter = 0;
    for (const string& file_name : files_list) {
        cout << "[" << (download_counter + 1) << "/" << files_list.size() 
             << "] Scaricando: " << file_name << "..." << endl;
        
        // Costruiamo il comando rclone copy usando il flag --include per prendere SOLO quel file specifico
        // Nota: usiamo le virgolette per gestire correttamente spazi nei percorsi di Linux
        string rclone_cmd = "rclone copy gdrive: \"" + target_folder + "\" "
                            "--drive-root-folder-id \"" + drive_id + "\" "
                            "--include \"" + file_name + "\" -P";
        
        // Eseguiamo il comando nel terminale di Linux tramite ROOT
        int status = gSystem->Exec(rclone_cmd.c_str());
        
        if (status == 0) {
            download_counter++;
        } else {
            cerr << " -> Errore durante il download di: " << file_name << endl;
        }
    }
    
    // 4. Rimozione dei file identificatori (Pulizia)
    cout << "\n=== Pulizia file identificatori di sistema ===" << endl;
    
    // Rimuove i file "._*" tipici dei metadati Apple/Drive o file temporanei nascosti
    string clean_cmd1 = "rm -f \"" + target_folder + "\"/._*";
    gSystem->Exec(clean_cmd1.c_str());
    
    // Opzionale: se per "file identifier" intendi file di testo di indice/configurazione rimasti,
    // puoi decommentare la riga sotto modificando l'estensione o il pattern da cancellare:
    // string clean_cmd2 = "rm -f \"" + target_folder + "\"/id_*.txt";
    // gSystem->Exec(clean_cmd2.c_str());
    
    cout << "\nProcesso completato! Scaricati con successo " << download_counter << " file su " << files_list.size() << "." << endl;
}