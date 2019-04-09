#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
ifstream fin("input.txt");
ofstream fout("output.txt");
class vector
{
protected:
    int nr;
    double *valori;
public:
    vector();
    vector(int n);
    vector(int n, double *v);
    vector(const vector& v);
    ~vector();
    virtual vector operator=(const vector &v);
    virtual vector operator+(const vector &v);
    virtual vector operator-(const vector &v);
    virtual vector operator*(const vector &v);
    //capacitate
    int dimensiune(); //intoarce numarul de valori
    void redimensionare(int n); //reduce/creste numarul elementelor de la nr la n
    bool neinitializat(); //verifica daca vectorul este neinitializat
    //element access
    double element(int poz); //intoarce elementul de pe pozitia poz
    double prim(); //intoarce primul element
    double ultim(); //intoarce ultimul element
    double* data(); //intoarce pointer catre primul element
    //modificatori
    void adaugaFinal(double x);//adauga element la sfarsit
    void stergeFinal();//sterge elementul de la sfarsit
    void insereaza(int poz,double x);//insereaza elementul x pe pozitia poz
    void stergePoz(int poz);
    void interschimba(vector& v);
    void sterge();
    //operatori relationali
    bool operator==(const vector &v);
    bool operator!=(const vector &v);
    bool operator<(const vector &v);
    bool operator>(const vector &v);
    bool operator>=(const vector &v);
    bool operator<=(const vector &v);
    //friends
    friend istream& operator>>(istream &i, vector& v);
    friend ostream& operator<<(ostream &o,const vector& v);
};

vector::vector()
{
    nr = 0;
    valori = new double;
    *valori = 0;
}
vector::vector(int n)
{
   nr = n;
    if(nr == -1 || nr == 0)
    {
        valori = new double;
        *valori = 0;
    }
    else
   {
        valori = new double[n];
        for(int i = 0; i < n; i++)
        valori[i] = 0;
   }
}
vector::vector(int n, double *v)
{
    nr = n;
    valori = new double[n];
    for(int i = 0; i < n; i++)
        valori[i] = v[i];
}
vector::vector(const vector& v)
{
    nr = v.nr;
    valori = new double[nr];
    for(int i = 0; i < nr; i++)
        valori[i] = v.valori[i];
}
vector::~vector()
{
    delete []valori;
}
vector vector::operator=(const vector& v)
{
    //operator de atribuire
    if(nr == 0) //polinomul nu a fost initializat decat in constructor
    {
        nr = v.nr;
        delete valori;
        valori = new double[nr];
    }
    if(nr != v.nr) //dimensiunea polinomului trebuie modificata
    {
        nr = v.nr;
        delete []valori;
        valori = new double[nr];
    }
    for(int i = 0; i < nr; i++)
            valori[i] = v.valori[i];
    return *this;
}
vector vector::operator+(const vector &v)
{
    vector s(fmax(v.nr,nr)); //creez un vector de lungime maxima
    int i = 0;
    while(i < v.nr && i < nr)
    {
        s.valori[i] = valori[i] + v.valori[i];
        i++;
    }
    while(i < v.nr)
    {
        s.valori[i] = v.valori[i];
        i++;
    }
    while(i < nr)
    {
        s.valori[i] = valori[i];
        i++;
    }
    return s;
}
vector vector::operator-(const vector &v)
{
    vector d(fmax(v.nr,nr)); //creez un vector de lungime maxima
    int i = 0;
    while(i < v.nr && i < nr)
    {
        d.valori[i] = valori[i] - v.valori[i];
        i++;
    }
    while(i < v.nr)
    {
        d.valori[i] = -v.valori[i];
        i++;
    }
    while(i < nr)
    {
        d.valori[i] = valori[i];
        i++;
    }
    return d;
}
vector vector::operator*(const vector &v)
{
    vector p(fmax(v.nr,nr)); //creez un vector de lungime maxima
    int i = 0;
    while(i < v.nr && i < nr) // cat timp avem 2, inmultesc, am 0 pana la capat
    {
        p.valori[i] = valori[i] * v.valori[i];
        i++;
    }
    return p;
}
int vector::dimensiune()
{
    return nr;
}
void vector::redimensionare(int n)
{
    vector cpy = *this;
    delete []this->valori;
    nr = n;
    valori = new double[nr];
    for(int i = 0; i < fmin(nr,cpy.nr); i++)
        valori[i] = cpy.valori[i];
    return;
}
bool vector::neinitializat()
{
    if(nr == 0)
        return true;
    return false;
}
double vector::element(int poz)
{
    if(poz < 0 || poz >= nr)
        return 0;
    return valori[poz];
}
double vector::prim()
{
    if(nr == 0)
        return 0;
    return valori[0];//exista pentru orice vector(din constructor)
}
double vector::ultim()
{
    if(nr==0) //vector gol
        return 0;
    return valori[nr-1];
}
double* vector::data()
{
    return valori;
}
void vector::adaugaFinal(double x)
{
    redimensionare(nr+1);
    valori[nr-1] = x;
    return;
}
void vector::stergeFinal()
{
    redimensionare(nr-1);
}
void vector::insereaza(int poz,double x)
{
    redimensionare(nr+1);
    for(int i = nr-1; i >= poz+1; i--)
        valori[i] = valori[i-1];
    valori[poz] = x;
    return;
}
void vector::stergePoz(int poz)
{
    for(int i = poz; i < nr-1; i++)
        valori[i] = valori[i+1];
    redimensionare(nr-1);
    return;
}
void vector::interschimba(vector& v)
{
    vector cpy = *this;
    *this = v;
    v = cpy;
    return;
}
void vector::sterge()
{
    vector cpy;
    *this = cpy;
    return;
}
bool vector::operator==(const vector &v)
{
    if(nr != v.nr)
        return false;
    for(int i = 0; i < nr; i++)
        if(valori[i] != v.valori[i])
        return false;
    return true;
}
bool vector::operator!=(const vector &v)
{
    if(*this == v)
        return false;
    return true;
}
bool vector::operator<(const vector &v)
{
    for(int i = 0; i < fmin(nr,v.nr); i++)
        if(valori[i] > v.valori[i])
        return false;
    return true;
}
bool vector::operator>(const vector &v)
{
    if(*this < v && *this != v)
        return true;
    return false;
}
bool vector::operator>=(const vector &v)
{
    if(*this > v || *this == v)
        return true;
    return false;
}
bool vector::operator<=(const vector &v)
{
    if(*this < v || *this == v)
        return true;
    return false;
}
istream& operator>>(istream &i, vector& v)
{
    i>>v.nr;
    delete []v.valori;
    v.valori = new double[v.nr];
    for(int k = 0; k < v.nr; k++)
        i>>v.valori[k];
    return i;
}
ostream& operator<<(ostream &o,const vector& v)
{
    if( v.nr == 0)
    {
        o<<"( )";
        return o;
    }
    o<<"(";
    if(v.nr >= 1)
    for(int k = 0; k < v.nr -1; k++)
    o<<v.valori[k]<<", ";
    o<<v.valori[v.nr-1];
    o<<")";
    return o;
}
void verifareAritmetice(vector *v1, int n)
{
    vector *v = new vector[n];
    for(int i = 0; i < n; i++)
        v[i] = v1[i];
    fout<<endl<<"Verificarea operatiilor aritmetice:"<<endl;
    int j = 0;
    for(int i = 0; i < n; i++)
    {
        fout<<v[i]<<" + "<<v[j]<<" = "<<v[i]+v[j]<<endl;
        fout<<v[i]<<" - "<<v[j]<<" = "<<v[i]-v[j]<<endl;
        fout<<v[i]<<" * "<<v[j]<<" = "<<v[i]*v[j]<<endl;
    }
    delete []v;
    return;
}
void verificareDim(vector *v1, int n)
{
    vector *v = new vector[n];
    for(int i = 0; i < n; i++)
        v[i] = v1[i];
    fout<<endl<<"Verificarea functiilor pentru dimensiune:"<<endl;
    for(int i = 0; i < n; i++)
    fout<<"Dimensiunea vectorului "<<v[i]<<" este:"<<v[i].dimensiune()<<"."<<endl;
    int j = 0;
    int dim1 = v[j].dimensiune()+1;
    int dim2 = v[j].dimensiune()-1;
    int dim3 = 0;
    fout<<"Redimensionez vectorul "<<v[j]<<" de dimensiune "<<v[j].dimensiune()<<" la dimensiunea "<<dim1<<"."<<endl;
    v[j].redimensionare(dim1);
    fout<<"\tVectorul rezultat: "<<v[j]<<endl;
    fout<<"Redimensionez vectorul "<<v[j]<<" de dimensiune "<<v[j].dimensiune()<<" la dimensiunea "<<dim2<<"."<<endl;
    v[j].redimensionare(dim2);
    fout<<"\tVectorul rezultat: "<<v[j]<<endl;
    fout<<"Redimensionez vectorul "<<v[j]<<" de dimensiune "<<v[j].dimensiune()<<" la dimensiunea "<<dim3<<"."<<endl;
    v[j].redimensionare(dim3);
    fout<<"\tVectorul rezultat: "<<v[j]<<endl;
    if(v[j].neinitializat() == true)
        fout<<"Vectorul este gol."<<endl;
    else
        fout<<"Vectorul nu este gol"<<endl;
    delete []v;
    return;
}
void verificareAcces(vector *v1, int n)
{
    vector *v = new vector[n];
    for(int i = 0; i < n; i++)
        v[i] = v1[i];
    fout<<endl<<"Verificarea functiilor de acces:"<<endl;
    int j = 0;
    int poz = v[j].dimensiune() - 1;
    fout<<"In vectorul "<<v[j]<<" avem pe pozitia "<<poz<<" elementul: "<<v[j].element(poz)<<"."<<endl;
    for(int i = 0; i < n; i++)
        fout<<"In vectorul "<<v[i]<<" primul element este "<<v[i].prim()<<", iar ultimul element este "<<v[i].ultim()<<"."<<endl;
    double val = 3.14;
    fout<<"In vectorul "<<v[j]<<" modific primul element cu valoarea "<<val<<"."<<endl;
    double *p = v[j].data();
    *p = val;
    fout<<"\tVectorul rezultat: "<<v[j]<<"."<<endl;
    delete []v;
    return;
}
void verificareModificatori(vector *v1, int n)
{
    vector *v = new vector[n];
    for(int i = 0; i < n; i++)
        v[i] = v1[i];
    fout<<endl<<"Verificarea modificatorilor:"<<endl;
    for(int i = 0; i < n; i++)
    {
            fout<<"In vectorul "<<v[i]<<" adaug la final valoarea "<<i*10<<"."<<endl;
            v[i].adaugaFinal(i*10);
            fout<<"\tVectorul rezultat: "<<v[i]<<"."<<endl;
            v[i].stergeFinal();
            fout<<"\tSterg valoarea adaugata si vectorul este: "<<v[i]<<"."<<endl;
    }
    if(n == 1)
        return;
    int j = 1;
    double val = 3.14;
    fout<<"Inserez in vectorul "<<v[j]<<" inainte de fiecare pozitie valoarea "<<val<<"."<<endl;
    int i = 0;
    int total = v[j].dimensiune();
    while(i < total)
    {
        v[j].insereaza(2*i,val);
        fout<<"\tLa pasul "<<i+1<<": "<<v[j]<<endl;
        i++;
    }
    fout<<"\tSterg elemenul de pe pozitia 0. "<<endl;
    v[j].stergePoz(0);
    fout<<"\t\t"<<v[j]<<endl;
    fout<<"\tSterg elemenul de pe pozitia "<<v[j].dimensiune()-2<<"."<<endl;
    v[j].stergePoz(v[j].dimensiune()-2);
    fout<<"\t\t"<<v[j]<<endl;
    fout<<"Interschimb vectorul v1:"<<v[0]<<" cu vectorul v2: "<<v[1]<<"."<<endl;
    v[0].interschimba(v[1]);
    fout<<"\tv1: "<<v[0]<<endl<<"\tv2: "<<v[1]<<endl;
    fout<<"Sterg vectorul v: "<<v[0]<<"."<<endl;
    v[0].sterge();
    fout<<"\tv: "<<v[0];
    delete []v;
}

class polinom:public vector
{
    int grad;
    // a mostenit vectorul coef si nr propriu-zis

public:
    polinom();
    polinom(int n);
    polinom(int n, double *v);
    polinom(const polinom& p);
    ~polinom();
    polinom operator=(const polinom& p);
    polinom operator+(const polinom& p);
    polinom operator-(const polinom& p);
    polinom operator*(const polinom& p);
    bool operator==(const polinom& p);
    double calcul_valoare(double x);
    friend istream& operator>>(istream& i,polinom& p);
    friend ostream& operator<<(ostream& o,const polinom& p);
    friend class pereche;
};

polinom::polinom()
{
    grad = 0;
    //este apelat automat constructorul vector
}
polinom::polinom(int n):vector(n+1)
{
    grad = n;
}
polinom::polinom(int n, double *v):vector(n+1,v)
{
    grad = n;
}
polinom::polinom(const polinom& p)
{
    grad = p.grad;
    nr = p.nr;
    for(int i = 0; i < nr; i++)
        valori[i] = p.valori[i];
}
polinom::~polinom()
{
    grad = 0;
}
polinom polinom::operator=(const polinom& p)
{
    grad = p.grad;
    if(nr == 0) //polinomul nu a fost initializat decat in constructor
    {
        nr = p.nr;
        delete valori;
        valori = new double[nr];
    }
    if(nr != p.nr) //dimensiunea polinomului trebuie modificata
    {
        nr = p.nr;
        delete []valori;
        valori = new double[nr];
    }
     valori = p.valori;
    return *this;
}
polinom polinom::operator+(const polinom &p)
{
    int grad_S; //gradul polinomului suma este maximul dintre gradele celuilalt
    if(grad > p.grad)
        grad_S = grad;
    else
        grad_S = p.grad;
    //creez un polinom de grad maxim cu coeficienti initial nuli
    int i = 0;
    polinom S(grad_S);
    while(i <= grad && i <= p.grad) //adun coeficient cu coeficient pana ajung la finalul polinomului de grad minim
    {
        S.valori[i] = valori[i] + p.valori[i];
        i++;
    }
    while(i <= grad) //daca primul polinom e de grad mai mare decat al doilea
    {
        S.valori[i] = valori[i];
        i++;
    }
     while(i <= p.grad) //daca al doilea polinom e de grad mai mare decat primul
    {
        S.valori[i] = p.valori[i];
        i++;
    }
    return S;
}
polinom polinom::operator-(const polinom &p)
{
    int grad_D; //gradul polinomului diferenta este initial maximul dintre gradele celuilalt
    if(grad == p.grad)
    {
        grad_D = grad;
        while(valori[grad_D] == p.valori[grad_D] && grad_D >=0) //daca se reduc termenii de grad maxim, scade gradul dif
        grad_D--;
        if(grad_D == -1) //daca gradul este -1, cele doua polinoame sunt identice, deci diferenta este polinomul nul
            {
                polinom D(-1);
                return D;
            }
    }
    else //daca cele doua polinoame nu au acelasi grad, gradul dif este  maximul dintre gradele celuilalt
    {
        if(grad > p.grad)
        grad_D = grad;
    else
        grad_D = p.grad;
    }
    polinom D(grad_D); //creez polinomul de grad maxim cu coeficienti initial nuli
    int i = 0;
    while(i <= grad && i <= p.grad && i <= grad_D) //scad coeficient cu coeficient pana ajung la finalul polinomului de grad minim, sau daca sunt egale, pana la gradul la care se reduc (grad_D)
    {
        D.valori[i] = valori[i] - p.valori[i];
        i++;
    }
    while(i <= grad && i <= grad_D) //daca primul polinom e de grad mai mare decat al doilea
    {
        D.valori[i] = valori[i];
        i++;
    }
     while(i <= p.grad && i <= grad_D) //daca al doilea polinom e de grad mai mare decat al doilea
    {
        D.valori[i] = -p.valori[i];
        i++;
    }
    return D;
}
polinom polinom::operator*(const polinom &p)
{
    int grad_P = grad + p.grad; //gradul polinomului produs este suma gradelor celor doua polinoame
    polinom P(grad_P); //creez un polinom de grad corespunzator cu coeficienti initial nuli
    for(int i = 0; i <= grad; i++)
        for(int j = 0; j <= p.grad; j++)
        P.valori[i+j] += valori[i]*p.valori[j]; //adun pe pozitia corespunzatoare produsul fiecarui coeficient din primul polinom cu fiecare coeficient din al doilea
    return P;
}
bool polinom::operator==(const polinom& p)
{
    if(grad != p.grad)
        return false;
    return vector::operator==(p);
}
double polinom::calcul_valoare(double x)
{
    if(grad == -1) //polinomul nul
        return 0;
    if(grad == 0) //polinomul format din termenul liber
        return *valori;
    double S = 0;
    for(int i = 0; i <= grad; i++)
        S+= pow(x,i)*valori[i];
    return S;
}
istream& operator>>(istream& i, polinom& p)
{
    i>>p.grad;
    p.nr = p.grad+1;
    delete []p.valori;
    p.valori = new double[p.nr];
    for(int k = 0; k < p.nr; k++)
        i>>p.valori[k];
    return i;
}
ostream& operator<<(ostream& o,const polinom& p)
{
    if(p.grad == -1)
        o<<0<<" ";
    else
    {
        polinom nul = polinom(p.grad); //polinomul cu coef nuli de grad egal cu cel de afisat
        if(nul == p) // daca polinomul are toti coeficientii nuli
        {
            o<<0<<" ";
            return o;
        }
        int start = 0;
        while(p.valori[start] == 0)
            start++; //pana la k, avem coeficienti nuli, deci nu afisam nimic
        if(start == 0) //daca termenul liber este nenul
        o<<p.valori[start]<<" ";
        else
        o<<p.valori[start]<<"x^"<<start<<" ";// primul monom nenul
    for(int k = start+1 ; k <= p.grad; k++)
        if(p.valori[k] == 1) //daca un coeficient este 1, il omit
        o<<"+"<<"x^"<<k<<" ";
        else
        if(p.valori[k] == -1) //daca un coeficient este -1, las doar semnul -
        o<<"-x^"<<k<<" ";
        else
        if(p.valori[k] > 0) //daca un coeficient este pozitiv, pun semnul + inainte
        o<<" +"<<p.valori[k]<<"x^"<<k<<" ";
        else
        if(p.valori[k] < 0) //daca un coeficient este negativ, semnul - apare oricum
        o<<p.valori[k]<<"x^"<<k<<" ";
    }
    return o;
}
class pereche
{
    double r;
    polinom p;
public:
    pereche();
    pereche(const polinom &p1, double d);
    pereche(const pereche &p1);
    bool verificare_radacina();
    bool operator==(const pereche& p);
    friend istream& operator>>(istream& i, pereche& pr);
    friend ostream& operator<<(ostream& i, const pereche& pr);
};
pereche::pereche()
{
    r = 0;
    polinom p(-1);
}
pereche::pereche(const polinom &p1, double d)
{
    r = d;
    p = p1;
}
pereche::pereche(const pereche &p1)
{
    r = p1.r;
    p = p1.p;
}
bool pereche::operator==(const pereche& pr)
{
  if(r != pr.r)
    return false;
  if(p == pr.r)
    return true;
  return false;
}

bool pereche::verificare_radacina()
{
    if(p.calcul_valoare(r) == 0) //daca polinomul p ia valoarea 0 pentru valoarea reala r, atunci e radacina
        return true;
    return false;
}

istream& operator>>(istream& i, pereche& pr)
{
    i>>pr.p;//citesc polinomul
    i>>pr.r;//citesc numarul real
    return i;
}

ostream& operator<<(ostream& i,const pereche& pr)
{
    i<<"polinomul:"<<pr.p<<"\t";
    i<<"numarul real:"<<pr.r;
    return i;
}

void verifOpAritm(polinom *p, int n)
{
    fout<<endl;
    for(int i = 0; i < n; i++)
        fout<<"Suma dintre polinomul "<<p[0]<<"si polinomul "<<p[i]<<": "<<p[0]+p[i]<<endl;
    fout<<endl;
    for(int i = 0; i < n; i++)
        fout<<"Diferenta dintre polinomul "<<p[0]<<"si polinomul "<<p[i]<<": "<<p[0]-p[i]<<"."<<endl;
    fout<<endl;
    for(int i = 0; i < n; i++)
        fout<<"Produsul dintre polinomul "<<p[0]<<"si polinomul "<<p[i]<<": "<<p[0]*p[i]<<"."<<endl;
    fout<<endl;
    for(int i = 0; i < n; i++)
        if(p[0] == p[i])
        fout<<"Polinomul "<<p[0]<<" este egal cu poliomul "<<p[i]<<"."<<endl;
    else
        fout<<"Polinomul "<<p[0]<<"nu este egal cu poliomul "<<p[i]<<"."<<endl;
    return;

}
int main()
{
    /*
    int n;
    fin>>n;
    vector *v;
    v = new vector[n];
    for(int i = 0; i < n; i++)
        fin>>v[i];
    for(int i = 0; i < n; i++)
        fout<<"Vectorul "<<i<<":  "<<v[i]<<endl;
    verifareAritmetice(v,n);
    verificareDim(v,n);
    verificareAcces(v,n);
    verificareModificatori(v,n);
    delete []v;
    */
    int n;
    fin>>n;
    polinom *p;
    p = new polinom[n];
    for(int i = 0; i < n; i++)
        fin>>p[i];
    for(int i = 0; i < n; i++)
        fout<<"polinomul "<<i<<":  "<<p[i]<<endl;
    verifOpAritm(p,n);
    delete []p;
    fin.close();
    fout.close();
    return 0;
}
