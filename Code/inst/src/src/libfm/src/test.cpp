#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
using namespace std;

//vector<vector<string> > main() {
int main() {
	ifstream infile;
	string str;
	infile.open("./code_testing_data.txt");

	if (!infile.is_open()) {
		cerr << "Unable to open file datafile.txt";
		exit(1);   // call system to stop
	}

	int total_line = 0;
	// include first line
	while(getline(infile, str)) {
		total_line ++;
	}

	infile.clear();
	infile.seekg(0);

	vector<vector<string> > ret;
	vector<string> vstr;
	vector<string> nstr;
	string::size_type sz;
	int row, col;
	int pos1, pos2;
	string item;	
	// from second line
	getline(infile, str);

	row = 0;
	while(getline(infile, str)) {
		pos1 = str.find("\t")+1;
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
	ret.push_back(vstr);
	ret.push_back(nstr);
	cout << ret[0][100] << endl;
	cout << ret[1][10] << endl;

	//return vstr;
	return 0;
}
