#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
// #include <boost/algorithm/string.hpp>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning
	vector<int> CN;		//gives number of occurence of each condition
public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
		nvalues=n;
		values=vals;
	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	vector<int> get_CN()
	{
		return CN;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
	void set_CN(vector<int> new_CN)
	{
		CN.clear();
		CN = new_CN;
	}
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }



};


 // The whole network represted as a list of nodes
class network{

	list <Graph_Node> Pres_Graph;

public:
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}
    
    
	int netSize()
	{
		return Pres_Graph.size();
	}
    // get the index of node with a given name
    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return count;
            count++;
        }
        return -1;
    }
	// get the node at nth index
    list<Graph_Node>::iterator get_nth_node(int n)
    {
       list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(count==n)
                return listIt;
            count++;
        }
        return listIt; 
    }
    //get the iterator of a node with a given name
    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return listIt;
        }
    
            cout<<"node not found\n";
        return listIt;
    }

	void print_graph()
	{
		list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            cout<<(listIt->get_name())<<endl;
			vector<string> p = listIt->get_Parents();
			vector<string> values = listIt->get_values();
			vector<float> cpt = listIt->get_CPT();
			cout<<p.size()<<" ";
			cout<<values.size()<<" "<<cpt.size();
			cout<<".";
			cout<<endl;
        }
	}

	void disply(){
		list<Graph_Node>::iterator listIt;
		for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++){
			string name = listIt->get_name();
			cout<<name<<" :"<<endl;
			vector<string> parL = listIt->get_Parents();
			// for (int q = 0; q < parL.size(); q++){
			// 	cout<<parL[q]<<", ";
			// }
			if(parL.size() == 0) cout<<"No parent"<<endl;
			else cout<<endl;
			vector<float> icpt = listIt->get_CPT();
			for (int q = 0; q < icpt.size(); q++){
				cout<<icpt[q]<<", ";
			}
			cout<<endl<<endl;
		}
	}

	void initalize()
	{
		list<Graph_Node>::iterator listIt;
		for (listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
		{
			vector<float> cpt = listIt->get_CPT();
			vector<string> values = listIt->get_values();
			int l1 = values.size();
			int l = cpt.size();
			for (int i = 0;i<l;i++)
			{
				cpt[i] = 1./l1;
			}
			listIt->set_CPT(cpt);
		}
	}

	bool update(vector<string> input)
	{
		int index = 0;
		int a = 0;
		bool out = true;
		list<Graph_Node>::iterator listIt;
		for (listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
		{
			string val = input[index];
			vector<string> values = listIt->get_values();
			vector<string> parents = listIt->get_Parents();
			int index_cpt = 0;
			int k = 1;
			// finding index whose probability needs to be increased

			for(int i = 0;i<parents.size();i++){
				list<Graph_Node>::iterator listIt1 = search_node(parents[parents.size()-i-1]);
				int p_index = get_index(parents[parents.size()-i-1]);
				vector<string> p_values = listIt1->get_values();
				string p_val = input[p_index];
				for (int l = 0;l<p_values.size();l++)
				{
					if (p_values[l] == p_val)
					{
						index_cpt += l*k;
						break;
					}
				}
				k = k*p_values.size();
			}
			vector<int> cn = listIt->get_CN();
			vector<float> cpt = listIt->get_CPT();
			int sum = 1;
			for (int l = 0;l<values.size();l++){
				sum += cn[index_cpt+l*k];
				if (values[l] == val)
					cn[index_cpt+l*k] +=1;
			}
			for (int l = 0;l<values.size();l++){
					float temp = cpt[index_cpt+l*k];
					cpt[index_cpt+l*k] = (float(cn[index_cpt+l*k]))/(float(sum));
					if (abs(temp - cpt[index_cpt+l*k])>0.000001)
						out = false;
			}
			listIt->set_CN(cn);
			listIt->set_CPT(cpt);
			index++;
		}
		return out;
	}

	vector<string> cpt_vector(vector<string> Vdata){
		int varIndex = 0;
		for (int i = 0; i < Vdata.size(); i++){
			if(Vdata[i] == "?"){
				varIndex = i;
				break;
			}
		}

		list<Graph_Node>::iterator itr = get_nth_node(varIndex);

		vector<string> values = itr->get_values();
		vector<float> cpt = itr->get_CPT();
		vector<string> parents = itr->get_Parents();
		
		int index_cpt = 0;
		int k = 1;
		for(int i = 0;i<parents.size();i++){

			list<Graph_Node>::iterator listIt1 = search_node(parents[parents.size()-i-1]);
			int p_index = get_index(parents[parents.size()-i-1]);
			vector<string> p_values = listIt1->get_values();
			string p_val = Vdata[p_index];
			for (int l = 0;l<p_values.size();l++)
			{
				if (p_values[l] == p_val)
				{
					index_cpt += l*k;
					break;
				}
			}
			k = k*p_values.size();
		}

		vector<float> ret;
		for (int l = 0;l<values.size();l++){
			ret.push_back(cpt[index_cpt+l*k]);
		}

		float max = 0.0;
		int indx = 0;
		for (int y = 0; y < ret.size(); y++){
			if(max <= ret[y]){
				max = ret[y];
				indx = y;
			}
		}
		Vdata[varIndex] = values[indx];
		return Vdata;
	}

	// markov blanket sampling
	vector<string> fill_val(vector<string> input)
	{
		list<Graph_Node>::iterator listIt;
		int index = 0;
		for (listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
		{

			string val = input[index];

			if (val == "?")
			{
				vector<string> values = listIt->get_values();
				vector<string> parents = listIt->get_Parents();
				vector<int> children = listIt->get_children();
				vector<float> cpt = listIt->get_CPT();
				int index_cpt = 0;
				int k = 1;
	// parent probability part
				for(int i = 0;i<parents.size();i++){
					list<Graph_Node>::iterator listIt1 = search_node(parents[parents.size()-i-1]);
					int p_index = get_index(parents[parents.size()-i-1]);
					vector<string> p_values = listIt1->get_values();
					string p_val = input[p_index];
					for (int l = 0;l<p_values.size();l++)
					{
						if (p_values[l] == p_val)
						{
							index_cpt += l*k;
							break;
						}
					}
					k = k*p_values.size();
				}
				vector<float> prob;
				for (int i = 0;i<val.size();i++)
				{
					int k1 = 1;
					int index_cpt1 = 0;
					float pi = cpt[index_cpt+i*k];
					string val1 = values[i];
					for (int j = 0;j<children.size();j++)
					{
						list<Graph_Node>::iterator listIt1 = get_nth_node(children[j]);
						vector<string> c_values = listIt1->get_values();
						vector<string> c_parents = listIt1->get_Parents();
						vector<float> c_cpt = listIt1->get_CPT();
				// child probability consideration part
						for(int i1 = 0;i1<c_parents.size();i1++){
							list<Graph_Node>::iterator listIt2 = search_node(parents[parents.size()-i1-1]);
							int p_index = get_index(c_parents[parents.size()-i1-1]);
							vector<string> p_values = listIt2->get_values();
							string p_val = input[p_index];
							for (int l = 0;l<p_values.size();l++)
							{
								if (p_values[l] == p_val)
								{
									index_cpt1 += l*k1;
									break;
								}
								else if (p_val=="?" && p_values[l] == val1)
								{
									index_cpt1 += l*k1;
									break;
								}
								
							}
							k1 = k1*p_values.size();
						}

						for (int i1 = 0;i1<c_values.size();i1++)
						{
							if (c_values[i1] == input[children[j]])
								index_cpt1 +=i1*k1;
						}
						pi *=c_cpt[index_cpt1];
					}
					prob.push_back(pi);
				}
				string out = values[0];
				float p_max = prob[0];
				// choosing best value
				for (int m = 0;m<prob.size();m++)
				{
					if (p_max<prob[m])
					{
						out = values[m];
						p_max = prob[m];
					}
				}
				input[index] = out;
			}
			index++;
		}
		return input;
	}

	void write_network(string filename){
		ifstream myfile(filename); 
		ofstream fileoo("solved_alarm.bif");
		string temp,line;
		
		while (! myfile.eof() ){
			stringstream ss;
			getline (myfile,line);
			ss.str(line);
			ss>>temp;
			if(! myfile.eof()) fileoo << line << endl;
			else {
				fileoo << line;
				break;
				}
			if(temp.compare("probability")==0){        
				ss>>temp;
				ss>>temp;
				// cout<< <<" ";
				list<Graph_Node>::iterator listIt = search_node(temp.substr(1,temp.length()-2));
				vector<float> cptVec = listIt->get_CPT();
				
				getline (myfile,line);
				stringstream ss2;
				ss2.str(line);
				ss2>> temp;
				string mainStr = "\t"+temp + " ";
				for (int i = 0; i < cptVec.size(); i++){
					stringstream ssj;
					if(cptVec[i] <= 0.00001) ssj << 0.00001;
					else ssj << cptVec[i];
					string x = ssj.str();
					mainStr += x+" ";
				}
				mainStr += ";";
				fileoo << mainStr << endl;
			}
		}
		myfile.close();
	}

};

network read_network(string filename){
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile(filename); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
					name = name.substr(1,name.length()-2);
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					
     					ss2>>temp;
     					
     					
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
						temp = temp.substr(1,temp.length()-2);
     					values.push_back(temp);
     					ss2>>temp;
    				}
					// cout<<name<<" ";
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);

     				
     		}
     		else if(temp.compare("probability")==0)
     		{       
				 	
     				ss>>temp;
     				ss>>temp;
     				temp = temp.substr(1,temp.length()-2);
                    list<Graph_Node>::iterator listIt;
                    list<Graph_Node>::iterator listIt1;
     				listIt=Alarm.search_node(temp);
                    int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
						temp = temp.substr(1,temp.length()-2);
                        listIt1=Alarm.search_node(temp);
						
                        listIt1->add_child(index);
     					values.push_back(temp);
     					
     					ss>>temp;

    				}
                    listIt->set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
					vector<int> curr_CN;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
     					curr_CPT.push_back(atof(temp.c_str()));
						curr_CN.push_back(0);	
     					ss2>>temp;
                       
    				}
                    
                    listIt->set_CPT(curr_CPT);
					listIt->set_CN(curr_CN);
     		}
            else
            {
                
            }		
    	}
    	
    	if(find==1)
    	myfile.close();
  	}
  	
  	return Alarm;
}

bool missing(vector<string> input)
{
	bool out = false;
	for (int i =0;i<input.size();i++)
		if (input[i] == "?")
			out = true;
	return out;
}
std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

int main(int argc,char* argv[])
{
	// argv[1] = "alarm.bif"
	// argv[2] = "records.dat"
	network Alarm;
	Alarm=read_network(argv[1]);
	Alarm.initalize();
	ifstream infile;
    infile.open(argv[2]);
	if (!infile.is_open()){
		cout<<".dat file not found"<<endl;
		return 0;
	}
	string X;
	vector< vector<string> > records;
    while(getline(infile,X)){
		vector<string> splitString;
    	// boost::split(splitString, X, boost::is_any_of(" "));
		splitString = split(X,' ');
		for (int i = 0; i < splitString.size(); i++){
			splitString[i] = splitString[i].substr(1,splitString[i].length()-2);
			// cout<<splitString[i]<<"/";
		}
		vector<string> expected;
		if (missing(splitString))
			expected = Alarm.cpt_vector(splitString);
		else 
			expected = splitString;
		bool a = true;
		while(a)
		{
			a = Alarm.update(expected);
		}
		records.push_back(splitString);
		// break;
    }
	// vector<string> out = Alarm.cpt_vector(records[0]);
	// for (int i = 0;i<out.size();i++)
	// 	cout<<out[i]<<" ";
	// 	cout<<""<<endl;
	// Alarm.disply();
	Alarm.write_network(argv[1]);
}