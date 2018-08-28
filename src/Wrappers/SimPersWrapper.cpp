/*
(c) 2016 Jyamiti Group, CSE, OSU
*/

#include "SimPersWrapper.hpp"

extern std::vector<int> complexSizes;
extern std::vector<int> accumulativeSizes;
extern int accumulativeSize;
extern float fThreshold;
extern vector<float> vecFiltrationScale;
extern SimplicialTree<bool> domain_complex;
//extern SimplicialTree<bool> range_complex;
//
//timer
extern clock_t start, timer1, timer2;
extern double dFuncTimeSum;
extern double dInsertTime;
extern double dCollapseTime;
extern int max_dimension;

void ComputingPersistenceForSimplicialMap(const char* file_name_of_domain_complex,
    bool is_domain_complex_with_annotation,
    const char* file_name_of_range_complex,
    const char* file_name_of_simplicial_map,
    bool is_save_range_complex_with_annotation = false,
    const char* new_range_complex_file_name = NULL);

void ComputingPersistenceForSimplicialMapElementary(const char* file_name_of_domain_complex,
    bool is_domain_complex_with_annotation,
    vector<string>& vecElemOpers,
    bool is_save_range_complex_with_annotation = false,
    const char* new_range_complex_file_name = NULL);

/***********************************************/

int SimpersWrapper::Run(vector<Operation *> &collapses,
                         //vector<Barcode> &output,
                         string fp) 

{
    start = clock();
    cout<<"Running Simpers"<<endl;
    std::string input_domain_complex_file_name;
    std::string input_range_complex_file_folder;
    std::string input_simplicial_map_file;      //file name for elementary mode and folder name for genaral mode
    std::string output_range_complex_with_annotation_file_name;
    std::string output_persistence_file_name;
    bool is_input_domain_complex_with_annotation = false;
    bool is_output_range_complex_with_annotation = false;
    bool is_elementary = true;
    //indicate if persistence barcode is ordered according to death time
    bool bDeathTimeOrder = true;
    int numFiltration;
    string sDelimiter("#");
    vecFiltrationScale.push_back(0);
    float fMaxScale = 0.0;
    input_domain_complex_file_name = fp+"_iDC";
    input_range_complex_file_folder = fp;
    input_simplicial_map_file = fp+"_collapses";
    output_range_complex_with_annotation_file_name = fp + "_range";
    output_persistence_file_name = fp + "_pers";
    is_input_domain_complex_with_annotation = false;
    is_output_range_complex_with_annotation = true;
    is_elementary = true;
    fThreshold = 0;
    numFiltration = 4;
    max_dimension = 2;
        
    
    ifstream ifs_m;
    ifs_m.open(input_simplicial_map_file);
    if (!ifs_m.is_open())
    {
        std::cout << "can't open file!" << endl;
        return 0;
    }
    filtration_step = 0;
    int stepsize = 0;
    ofstream myfile;
    myfile.open ("fromRandom.txt",ios::out);
    auto t = collapses[stepsize];
    while( stepsize<collapses.size() ) //while (!(ifs_m.eof()))
    {
        //t= collapses[stepsize];
        vector<string> vecOpers;
        string sOper;
        char sLine[256];
        bool bValid = false;
        while(t = collapses[stepsize])//while (ifs_m.getline(sLine, 256))
        {
            
            strcpy(sLine,(t->PrintString()).c_str());
            myfile<<sLine;
            if (sLine[0] == '#')
            {
                bValid = true;
                stringstream ss;
                ss.str(sLine);
                string sSharp;
                float fScale;
                ss >> sSharp;
                ss >> fScale;
                vecFiltrationScale.push_back(fScale);
                break;
            }
            sOper = sLine;
            vecOpers.push_back(sOper);
            stepsize++;
        }
        if (bValid)
        {
            //timer
            timer2 = clock();
            ComputingPersistenceForSimplicialMapElementary(input_domain_complex_file_name.c_str(),
                is_input_domain_complex_with_annotation,
                vecOpers, false,
                output_range_complex_with_annotation_file_name.c_str());
            dFuncTimeSum += (clock() - timer2);
        }
        stepsize++;
    }
    //cout<<"Exited loop check for the from Random file";
    //getchar();
    ifs_m.close();
    ifs_m.clear();
    myfile.close();
    myfile.clear();
    if (is_output_range_complex_with_annotation)
        domain_complex.WriteComplexWithAnnotation(output_range_complex_with_annotation_file_name.c_str());

    ofstream ofsTime("Timing.txt");
    cout << "Simpers Total time without I/O: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
    ofsTime << "Simpers Total time without I/O: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
    //output persistence barcode
    ofstream ofs_dgm(output_persistence_file_name);
    for (int i = 0; i < persistences.size(); ++i)
    {
        std::vector<std::pair<int, int>> vecBarcodes;
        for (std::unordered_map<int, pair<int, int>>::iterator it = persistences[i].begin(); it != persistences[i].end(); ++it)
        {
            if (bDeathTimeOrder)
                vecBarcodes.push_back(std::make_pair((it->second).first, (it->second).second));
            else
            {
                ofs_dgm << i << " " << vecFiltrationScale[(it->second).first] << " ";
                if ((it->second).second == -1)
                    ofs_dgm << "inf" << endl;
                else
                    ofs_dgm << vecFiltrationScale[(it->second).second] << endl;
            }
        }
        if (bDeathTimeOrder)
        {
            std::stable_sort(vecBarcodes.begin(), vecBarcodes.end(), barcodeCompare);
            for (int k = 0; k < vecBarcodes.size(); ++k)
            {
                ofs_dgm << i << " " << vecFiltrationScale[vecBarcodes[k].first] << " ";
                if (vecBarcodes[k].second == -1)
                    ofs_dgm << "inf" << endl;
                else
                    ofs_dgm << vecFiltrationScale[vecBarcodes[k].second] << endl;
            }
        }
    }
    ofs_dgm.close();
    ofs_dgm.clear();
    
    //find the maximum scale
    for (int k = 0; k < vecFiltrationScale.size(); ++k)
    {
        if (vecFiltrationScale[k] > fMaxScale)
            fMaxScale = vecFiltrationScale[k];
    }

    //output timing result
    dFuncTimeSum = dFuncTimeSum / (double)CLOCKS_PER_SEC;
    
    ofsTime << "Simpers Total time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
    cout << "Simpers Total time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
    //ofsTime << "Simpers computation time: " << dFuncTimeSum << "s" << endl;
    //cout << "Simpers computation time: " << dFuncTimeSum << "s" << endl;

    ofstream ofs_size("size.txt");
    int maxCS = 0;
    ofs_size << "Scales : Complex Sizes" << endl;
    for (int i = 0; i < vecFiltrationScale.size(); ++i)
    {
        ofs_size << vecFiltrationScale[i] << " : " << complexSizes[i] << endl;
        if (complexSizes[i] > maxCS)
            maxCS = complexSizes[i];
    }
    ofs_size << "Cumulative Complex size: " << accumulativeSizes.back() << endl;
    ofs_size << endl;
    ofs_size.close();
    ofs_size.clear();

    ofsTime << "Max Complex Size: " << maxCS << endl;
    cout << "Max Complex Size: " << maxCS << endl;
    ofsTime << "Cumulative Complex Size: " << accumulativeSizes.back() << endl;
    cout << "Cumulative Complex Size: " << accumulativeSizes.back() << endl;
    ofsTime << "Max Scale: " << fMaxScale << endl;
    cout << "Max Scale: " << fMaxScale << endl;
    ofsTime.close();
    ofsTime.clear();

return 0;
}