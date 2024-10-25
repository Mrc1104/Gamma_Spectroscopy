#ifndef STRTOOLS_H
#define STRTOOLS_H

// parse string with vars seperated by a delim
vector<string> parse_str(string &str, char delim = ',')
{
	if(str.back() != delim)
		str += delim;
		
	std::vector<std::string> sfiles;
	char trash[755]; int ti = 0;
	memset(trash, 0, 255);
	for(int i = 0; i < str.size(); i++) {
		if( str[i] == delim ){
			sfiles.push_back(trash);
			// std::cout << sfiles.back() << std::endl;
			memset(trash, 0, 255);
			ti = 0;
		}
		else
			trash[ti++] = str[i];
	}
	return sfiles;
}

#endif
