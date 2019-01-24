#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
using namespace std;


class PrePostProcessor {
    public:
        PrePostProcessor() {};
        vector<vector<string> > preprocessor(string filepath);
        vector<vector<string> > preprocessor2(vector<string> data_table);
        vector<string> postprocessor(vector<double> vpred);
        std::vector<std::string> obs_lines;
    protected:
        int total_line;
};

vector<vector<string> > PrePostProcessor::preprocessor2(vector<string> data_table) {
	string str;

	total_line = data_table.size();
	// include first line

	vector<vector<string> > ret;
	vector<string> vstr;
	vector<string> nstr;
    vector<string> vmeta;
	string::size_type sz;
	int row, col;
	// int pos_first_tab, pos1, pos2;
  std::string::size_type pos_first_tab, pos1, pos2;
	string item;	
	// from second line
	str = data_table[0];
	obs_lines.push_back(str);

	row = 0;
  col = total_line;
	//while(getline(infile, str)) {
	for(int data_i=1; data_i<total_line; data_i++) {
	  str = data_table[data_i];
        obs_lines.push_back(str);
        pos_first_tab = str.find("\t");

		pos1 = pos_first_tab+1;
		col = total_line;
		while((pos2 = str.find("\t", pos1)) != string::npos) {
			item = str.substr(pos1, pos2-pos1);
			try {
				stof(item, &sz);
				item += " " + to_string(row) + ":1 " + to_string(col) + ":1";
				vstr.push_back(item);
			}
			catch (const invalid_argument& ia) {
				item = "0 " + to_string(row) + ":1 " + to_string(col) + ":1";
				nstr.push_back(item);
			}

			pos1 = pos2 + 1;
			col += 1;
		}
		item = str.substr(pos1, pos2-pos1);
		try {
			stof(item, &sz);
			item += " " + to_string(row) + ":1 " + to_string(col) + ":1";
			vstr.push_back(item);
		}
		catch (const invalid_argument& ia) {
			item = "0 " + to_string(row) + ":1 " + to_string(col) + ":1";
			nstr.push_back(item);
		}
		row += 1;
	}
    int i; 
    for(i=0; i<total_line; i++) {
        vmeta.push_back("0");
    }
    for(i=total_line; i<col; i++) {
        vmeta.push_back("1");
    }

	ret.push_back(vstr);
	ret.push_back(nstr);
    ret.push_back(vmeta);

    cout << total_line << endl;
    cout << col << endl;
	cout << ret[0][100] << endl;
	cout << ret[1][10] << endl;

	return ret;
	//return 0;
}


vector<vector<string> > PrePostProcessor::preprocessor(string filepath) {
  ifstream infile;
  string str;
  //infile.open("./code_testing_data.txt");
  
  infile.open(filepath.c_str());
  
  if (!infile.is_open()) {
    cerr << "Unable to open file datafile.txt";
    exit(1);   // call system to stop
  }
  
  total_line = 0;
  // include first line
  while(getline(infile, str)) {
    total_line ++;
  }
  
  infile.clear();
  infile.seekg(0);
  
  vector<vector<string> > ret;
  vector<string> vstr;
  vector<string> nstr;
  vector<string> vmeta;
  string::size_type sz;
  int row, col;
  std::string::size_type pos_first_tab, pos1, pos2;
  string item;	
  // from second line
  getline(infile, str);
  obs_lines.push_back(str);
  
  row = 0;
  col = total_line;
  while(getline(infile, str)) {
    obs_lines.push_back(str);
    pos_first_tab = str.find("\t");
    
    pos1 = pos_first_tab+1;
    col = total_line;
    while((pos2 = str.find("\t", pos1))!= string::npos) {
      item = str.substr(pos1, pos2-pos1);
      try {
        stof(item, &sz);
        item += " " + to_string(row) + ":1 " + to_string(col) + ":1";
        vstr.push_back(item);
      }
      catch (const invalid_argument& ia) {
        item = "0 " + to_string(row) + ":1 " + to_string(col) + ":1";
        nstr.push_back(item);
      }
      
      pos1 = pos2 + 1;
      col += 1;
    }
    item = str.substr(pos1, pos2-pos1);
    try {
      stof(item, &sz);
      item += " " + to_string(row) + ":1 " + to_string(col) + ":1";
      vstr.push_back(item);
    }
    catch (const invalid_argument& ia) {
      item = "0 " + to_string(row) + ":1 " + to_string(col) + ":1";
      nstr.push_back(item);
    }
    row += 1;
  }
  infile.close();
  int i; 
  for(i=0; i<total_line; i++) {
    vmeta.push_back("0");
  }
  for(i=total_line; i<col; i++) {
    vmeta.push_back("1");
  }
  
  ret.push_back(vstr);
  ret.push_back(nstr);
  ret.push_back(vmeta);
  
  cout << total_line << endl;
  cout << col << endl;
  cout << ret[0][100] << endl;
  cout << ret[1][10] << endl;
  
  return ret;
  //return 0;
}



vector<string> PrePostProcessor::postprocessor(vector<double> vpred) {
    
    vector<string> ret; 
	string::size_type sz;
	int row, col;
	std::string::size_type pos_first_tab, pos1, pos2;
	string item;	


	int pred_index = 0;
    // int i;
    std::string::size_type i;
    string header = obs_lines[0];
    ret.push_back(header);

    for(i=1; i<obs_lines.size(); i++) {
        pos_first_tab = obs_lines[i].find("\t");
        
        string line_pivot = obs_lines[i].substr(0, pos_first_tab);

		pos1 = pos_first_tab+1;
		col = total_line;
		while((pos2 = obs_lines[i].find("\t", pos1))!= string::npos) {
			item = obs_lines[i].substr(pos1, pos2-pos1);
			try {
				stof(item, &sz);
                line_pivot += string("\t") + item;
			}
			catch (const invalid_argument& ia) {
                line_pivot += string("\t") + to_string(vpred[pred_index]);
                pred_index++;
			}

			pos1 = pos2 + 1;
			col += 1;
		}
		item = obs_lines[i].substr(pos1, pos2-pos1);
		try {
			stof(item, &sz);
            line_pivot += string("\t") + item;
		}
		catch (const invalid_argument& ia) {
            line_pivot += string("\t") + to_string(vpred[pred_index]);
            pred_index++;
		}
		row += 1;
        ret.push_back(line_pivot);
	}

    return ret;
}
