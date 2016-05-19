#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <bitset>
#include <ios>
#include <iomanip>

using namespace std;

bitset<16> compress(fstream& in, fstream& out) {

	int i=0, h=1, sgsh=0, sty=0, n = 1;
	long t = 314, ty=0, b=0;
	vector<bitset<16>> x, xta, jxta;
	vector<bool> tnode,ynode;
	vector<bool> vnode, bnode;
	vector<unsigned char> uch,ych,tch,vch,bch;
	bitset<64> last;
	bitset<8> temp2, kbp, pbk;
	bitset<2> upg;
	vector<bitset<8>> page;
	vector<bitset<2>> upage;
	bitset<16> back;
	in.seekg(in.beg,in.end);
	long length = in.tellg();
	in.seekg(in.end,in.beg);
	long u, j = 0;

	if (length%2 != 0)
		n = 0;

	if (length < 50) {
		cout << "File assumed too small for compression..";
		exit(1);
	}
	while (h) {	
		// Clear byte (z)
		bitset<16> v, h3, h2, h1, c1, c2 = 0;
		int a=0, i=0, m=0;
		vector<bitset<16> > k, y;
		char l = '\0';

		// Create one pair of keys for each tier

		do {
			i = 0;
			do {
				t++;
				k.push_back(t);
				i++;
			} while (i < 2);

			h3.reset();
		
			// Cycle Keys and Create final byte
			for (int g = 0;g<2;g++) {
				// 'xta' will always have 'h3's in it. That's a non-debatable
				// This entire 'h' while loops till xta only has one left. (top of pyramid)
				if (j == 0) {
					in.get(l);
					y.push_back(static_cast<int>(l));
				}
				else {
					l = xta[0].to_ulong();
					if (xta.size() > 1)
						xta.erase(xta.begin()+1);
					else
						y.push_back(0);
					y.push_back(static_cast<long>(l));
				}
			}
			if (y.size() < 2)
				y.push_back(0);
			h1 = (y[0] ^ k[0]);
			h2 = (h1 ^ y[1]);
			h3 = (h2 ^ k[1]);
			jxta.push_back(h3);

			// y0 ^ k0 = h1
			// h1 ^ y1 = h2
			// h2 ^ k1 = h3
			// find solution for binary steps.
			// put in 1x1x1-bit binary for HEX bin tree output
			for (b = 1; b >= 0; b--) {
				back = y[b];
				if (b == 1) {
					v = (h3 ^ k[1] ^ k[0]);
					v = (h3 ^ v);
				} else {
					if (v == y[0]) {
						c2 = v;
						h2 = (h3 ^ k[1]);
						h1 = (c2 ^ h2);
						c1 = (h1 ^ k[0]);
						vnode.push_back(0);
						tnode.push_back(1);
						ty++;
						break;
					}
					else if (v == y[1]) {
						c1 = v;
						h1 = (c1 ^ k[0]);
						h2 = (h3 ^ k[1]);
						c2 = (h1 ^ h2);
						vnode.push_back(1);
						tnode.push_back(1);
						ty++;
						break;
					}
				}

				u = 0;

				for (int cnt = 0;cnt < 8;cnt++) {
					kbp[cnt] = back[cnt+8];
					pbk[cnt] = back[cnt];
				}

				if (pbk.to_string() != v.to_string().substr(0,7)) {
					for (int g=0;g<8;g++)
						temp2[g] = v[g];

					if (pbk.all()) {
						ty++;
						bnode.push_back(1);
						tnode.push_back(0);
						u = 1;
					}
					if (pbk.none()) {
						ty++;
						bnode.push_back(0);
						tnode.push_back(0);
						u = 1;
					}
				}
				else tnode.push_back(1);

				//cout << kbp.to_string() << endl << v.to_string().substr(32,63) << endl << endl;
				if (kbp.to_string() != v.to_string().substr(8,15)) {
					for (int g=8;g<16;g++)
						temp2[g] = v[g];
					if (kbp.all() == 1) {
						ty++;
						bnode.push_back(1);
						tnode.push_back(0);
						u = u + 2;
					}
					if (kbp.none() == 1) {
						ty++;
						bnode.push_back(0);
						tnode.push_back(0);
						u = u + 2;
					}
				}
				else tnode.push_back(1);

				if (u == 1)
					page.push_back(kbp);
				if (u == 2)
					page.push_back(pbk);
				if (u == 0) {
					page.push_back(kbp);
					page.push_back(pbk);
				}
				upage.push_back(u);
				u = 0;
					
				
			}
			// Clear keys, save byte
			y.clear();

			if (tnode.size() >= 64) {
				for (int g = 0;g<64;g++)
					last[g]=tnode[g];
				tch.push_back(static_cast<unsigned char>(last.to_ulong()));
				tnode.erase(tnode.begin(),tnode.begin()+64);
			}

			if (vnode.size() >= 64) {
				for (int g = 0;g<64;g++)
					last[g]=vnode[g];
				vch.push_back(static_cast<unsigned char>(last.to_ulong()));
				vnode.erase(vnode.begin(),vnode.begin()+64);
			}

			if (bnode.size() >= 64) {
				for (int g = 0;g<64;g++)
					last[g]=bnode[g];
				bch.push_back(static_cast<unsigned char>(last.to_ulong()));
				bnode.erase(bnode.begin(),bnode.begin()+64);
			}

			if (upage.size() >= 32) {
				b = 0;
				int op=0;
				upg = upage[0];
				while (op < upage.size()) {
					if (b%2 == 0) {
						op++;
						last[b%64] = upg[0];
						last[(b+1)%64] = upg[1];
						upg = upage[op];
					}
					if (b%64 == 0) {
						uch.push_back(static_cast<unsigned char>(last.to_ulong()));
						upage.erase(upage.begin(),upage.begin()+32);
					}
					b++;
				}
			}

		} while(in.peek() != EOF && xta.size() != 1);

		// Reset vars, move file position (see above)
		// swap out jxta/xta to set len of Keys Cycle
		n = 1;
		xta.clear();
		xta.swap(jxta);
		jxta.clear();
		j = xta.size();	

		if (j < 2)
			break;
	}

	for (int i=0;i<tch.size();i++)
		out << tch[i];

		for (int i=0;i<tnode.size();i++) {
			last[i%64] = tnode[i];
			if (i%64 == 0) {
				out << endl << static_cast<char>(last.to_ulong()) << flush;
				last.reset();
			}
		}
		tnode.clear();
		out << static_cast<char>(last.to_ulong()) << flush;
		last.reset();
		out << endl << "t" << endl << flush;

	for (int i=0;i<vch.size();i++)
		out << vch[i];

		for (int i=0;i<vnode.size();i++) {
			last[i%64] = vnode[i];
			if (i%64 == 0) {
				out << static_cast<char>(last.to_ulong()) << flush;
				last.reset();
			}
		}
		vnode.clear();
		out << static_cast<char>(last.to_ulong()) << flush;
		last.reset();
		out << endl << "v" << endl << flush;

	for (int i=0;i<bch.size();i++)
		out << bch[i];

		for (int i=0;i<bnode.size();i++) {
			last[i%64] = bnode[i];
			if (i%64 == 0) {
				out << static_cast<char>(last.to_ulong()) << flush;
				last.reset();
			}
		}
		bnode.clear();
		out << static_cast<char>(last.to_ulong()) << flush;
		last.reset();
		out << endl << "b" << endl << flush;

	int op = 0;
	for (int i=0;i<uch.size();i++)
		out << uch[i];

		if (upage.size() >= 1)
			upg = upage[0];
		while (upage.size() > op) {
			if (b%2 == 0) {
				op++;
				upg = upage[op];
			}
			if (b%64 == 0)
				out << static_cast<char>(last.to_ulong()) << flush;
			last[b%64] = upg[b%2];
			b++;
		}
		upage.clear();
		out << static_cast<char>(last.to_ulong()) << flush;
		out << endl << "u" << endl << flush;

	j=0;
	if (page.size() >= 1)
		temp2 = page[0];
	while (page.size() > j) {
		if (b%8 == 0) {
			j++;
			temp2 = page[j];
		}
		if (b%64 == 0)
			out << static_cast<char>(last.to_ulong()) << flush;
		last[b%64] = temp2[b%8];
		b++;
	}
	page.clear();
	out << static_cast<char>(last.to_ulong()) << flush;
	out << endl << "p" << endl << flush;

	out << n << endl << "n" << endl;

	cout << "Within: " << sty << " Without: " << sgsh << endl << ty;
	bitset<16> hey = xta[0];
	out.seekp(in.end,in.beg);
	out << static_cast<char>(t);
	return hey;
}

bool getVect(vector<bitset<64>> node, int pos) {
	static long y;
	bitset<64> n;
	n = node[node.size() - y - 1];
	return static_cast<bool>(n[pos]);
	if (pos == 0)
		y++;
}
			

void decomp(fstream& in, fstream& out) {
	int j = 0;
	string s, en;
	bitset<64> vtemp, ttemp, ptemp, btemp;
	bitset<64> hashed;
	bitset<8> bk, kb;
	vector<bitset<16>> kbp, pbk;
	vector<bool> tn;
	vector<int> swtch;
	vector<bitset<64>> tnode, vnode, bnode, page, upage;
	bitset<16> v, h3, h2, h1, c1, c2 = 0;
	vector<bitset<16>> k, y;
	bitset<64> h, n;
	bitset<2> d;
	vector<bitset<2>> upf;
	bitset<8> cxr,cxl;
	vector<bitset<16>> byte;
	bool TRUE=1, FALSE=0;
	int c=0,m=0,i=0;
	in.seekg(in.beg,in.end);
	long length = in.tellg();
	in.seekg(in.end,in.beg);
	do {
		getline(in, s, '\n');
		getline(in, en, '\n');
		if (en == "t") {
			for (char i : s)
				tnode.push_back(static_cast<int>(i));
		}
		if (en == "v") {
			for (char i : s)
				vnode.push_back(static_cast<int>(i));
		}
		if (en == "b") {
			for (char i : s)
				bnode.push_back(static_cast<int>(i));
		}
		if (en == "p") {
			for (char i : s)
				page.push_back(static_cast<int>(i));
		}
		if (en == "u") {
			for (char i : s)
				upage.push_back(static_cast<int>(i));
		}
	} while (en != "u");
	

	while (upage.size() - c - 1 >= 0) {
		n = upage[upage.size() - c - 1];
		for (int g=63;g>=0;g--) {
			d[g%2] = n[g];
			if (g%2 == 0)
				upf.push_back(d);
		}
		c++;
	}
	c = 1;
	
	for (int t=0;t<upf.size();t++) {
		v = (h3 ^ k[0] ^ k[1]);
		v = (h3 ^ v);
		k.erase(k.begin(),k.begin()+2);
		if (c%2 == 0 && getVect(tnode,m) == TRUE) {
			if (getVect(vnode,vi) == FALSE) {
				m++;
				vi++;
				c2 = v;
				h2 = (h3 ^ k[1]);
				h1 = (c2 ^ h2);
				c1 = (h1 ^ k[0]);
				c++;
			}
			else {
				m++;
				vi++;
				c1 = v;
				h1 = (c1 ^ k[0]);
				h2 = (h3 ^ k[1]);
				c2 = (h1 ^ h2);
				c++;
			}
			t--;
			continue;
		}				
		if (upf[t].to_ulong() == 0 && getVect(tnode,m) == TRUE && getVect(tnode,m+1) == TRUE) {
			m += 2;
			byte.push_back(v);
		}
		else if (upf[t].to_ulong() == 0 && getVect(tnode,m) == TRUE && getVect(tnode,m+1) == FALSE) {
			m+=2;
			if (getVect(bnode,ni) == TRUE) {
				bk.set();
				ni++;
				bitset<16> von = bk.to_string() + v.to_string().substr(0,7);
			}
			else {
				kb.reset();
				ni++;
				bitset<16> von =  kb.reset().to_string() + v.to_string().substr(0,7);
			}
			byte.push_back(von);
		}
		else if (upf[t].to_ulong() == 0 && getVect(tnode,m) == FALSE && getVect(tnode,m+1) == TRUE) {
			m += 2;
			if (getVect(bnode,ni) == TRUE) {
				bk.set();
				bitset<16> von = v.to_string().substr(8,15) + bk.to_string();
				byte.push_back(von);
			}
			else if (getVect(bnode,ni) == FALSE) {
				kb.reset();
				bitset<16> von = v.to_string().substr(8,15) + kb.to_string();
				byte.push_back(von);
			}
		}
		else if (upf[t].to_ulong() == 0 && getVect(tnode,m) == FALSE && getVect(tnode,m+1) == FALSE) {
			m += 2;
			if (getVect(bnode,ni) == TRUE) {
				bk.set();
				bitset<16> von = v.to_string().substr(8,15) + bk.to_string();
				byte.push_back(von);
			}
			else if (getVect(bnode,ni) == FALSE) {
				kb.reset();
				bitset<16> von = v.to_string().substr(8,15) + kb.to_string();
				byte.push_back(von);
			}
		}
		else {
			if (upf[t].to_ulong() == 1) {
				if (getVect(bnode,ni) == FALSE)
					cxr.reset();
				else
					cxr.set();
				ni++;
			}
			if (upf[t].to_ulong() == 2) {
				if (getVect(bnode,ni) == FALSE)
					cxl.reset();
				else
					cxl.set();
				ni++;
			}
			if (upf[t].to_ulong() == 3) {
				if (getVect(bnode,ni) == FALSE) {
					cxl.reset();
					if (getVect(bnode,ni) == FALSE)
						cxr.reset();
					else
						cxr.set();
				}
				else {
					cxl.set();
					bitset<16> von = cxl.to_string() + v.to_string();
				}
				byte.push_back(von);
				ni++;
			}
		}
		c++;
	}
	page.clear();
}

int main(int c, char * argv[]) {

	const char * m = argv[1];
	const char * n = argv[2];

	fstream in { m, std::ios_base::in | std::ios_base::binary };
	if (! in)
		{ cout << "no input file\n"; return -1; }

	fstream out { n, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary };
	if (! out)
		{ cout << "no output file\n"; return -1; }

	out.close();

	out.open(n, std::ios_base::out | std::ios_base::app | std::ios_base::binary);
	if (! out)
		{ cout << "no output file\n"; return -1; }

	cout << "Comb Compressor - 2016\n";

	bitset<16> xt = compress(in,out);
	out.seekp(out.end,out.beg);
	out << xt << std::endl;

	cout << "Finished!";

	in.close();
	out.close();
	return 0;
}