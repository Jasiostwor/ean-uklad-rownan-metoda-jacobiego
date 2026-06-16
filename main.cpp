#include <wx/wx.h>
#include <wx/valnum.h>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include "Jacobi.h"

using namespace std;
using namespace interval_arithmetic;

class JacobiApp : public wxApp {
public:
    virtual bool OnInit();
};

class JacobiFrame : public wxFrame {
public:
    JacobiFrame(const wxString& title);

private:
    wxTextCtrl* nInput;
    wxTextCtrl* mitInput;
    wxTextCtrl* epsInput;
    
    wxTextCtrl* inputMemo;
    wxTextCtrl* outputMemo;
    wxRadioBox* modeRadio;
    wxChoice* templateChoice;
    wxButton* solveButton;
    wxButton* templateButton;

    void OnSolve(wxCommandEvent& event);
    void OnTemplate(wxCommandEvent& event);
    
    struct ParsedData {
        int n = 0;
        int mit = 100;
        long double eps = 1e-14;
        vector<pair<int, int>> a_indices;
        vector<string> a_values;
        vector<int> b_indices;
        vector<string> b_values;
        vector<int> x_indices;
        vector<string> x_values;
    };
    
    bool ParseInputText(ParsedData& data);
    string Trim(const string& str);
};

wxIMPLEMENT_APP(JacobiApp);

bool JacobiApp::OnInit() {
    Interval<long double>::Initialize();
    
    JacobiFrame* frame = new JacobiFrame("Rozwiązywanie układów równań - Metoda Jacobiego");
    frame->Show(true);
    return true;
}

JacobiFrame::JacobiFrame(const wxString& title)
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(900, 700)) {
    
    wxPanel* panel = new wxPanel(this, wxID_ANY);
    
    wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
    
    // Top Parameters Sizer
    wxBoxSizer* paramSizer = new wxBoxSizer(wxHORIZONTAL);
    
    paramSizer->Add(new wxStaticText(panel, wxID_ANY, "n:"), 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    nInput = new wxTextCtrl(panel, wxID_ANY, "4");
    paramSizer->Add(nInput, 0, wxALL, 5);
    
    paramSizer->Add(new wxStaticText(panel, wxID_ANY, "mit:"), 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    mitInput = new wxTextCtrl(panel, wxID_ANY, "100");
    paramSizer->Add(mitInput, 0, wxALL, 5);
    
    paramSizer->Add(new wxStaticText(panel, wxID_ANY, "eps:"), 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    epsInput = new wxTextCtrl(panel, wxID_ANY, "1e-14");
    paramSizer->Add(epsInput, 0, wxALL, 5);
    
    mainSizer->Add(paramSizer, 0, wxEXPAND | wxALL, 5);
    
    wxBoxSizer* topSizer = new wxBoxSizer(wxHORIZONTAL);
    
    wxBoxSizer* inputSizer = new wxBoxSizer(wxVERTICAL);
    inputSizer->Add(new wxStaticText(panel, wxID_ANY, "Dane wejściowe (a[i,j], b[i], x[i]):"), 0, wxALL, 5);
    inputMemo = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(500, 200), wxTE_MULTILINE);
    inputSizer->Add(inputMemo, 1, wxEXPAND | wxALL, 5);
    
    wxBoxSizer* rightSizer = new wxBoxSizer(wxVERTICAL);
    wxArrayString modes;
    modes.Add("Zwykła zmiennopozycyjna");
    modes.Add("Przedziałowa (dane rzecz.)");
    modes.Add("Przedziałowa (dane przedz.)");
    modeRadio = new wxRadioBox(panel, wxID_ANY, "Tryb obliczeń", wxDefaultPosition, wxDefaultSize, modes, 1, wxRA_SPECIFY_COLS);
    
    templateChoice = new wxChoice(panel, wxID_ANY);
    templateChoice->Append("Szablon domyślny (zera)");
    templateChoice->Append("Przykład a)");
    templateChoice->Append("Przykład b)");
    templateChoice->Append("Przykład c)");
    templateChoice->SetSelection(0);
    
    templateButton = new wxButton(panel, wxID_ANY, "Generuj szablon", wxDefaultPosition, wxSize(150, 30));
    templateButton->Bind(wxEVT_BUTTON, &JacobiFrame::OnTemplate, this);
    
    solveButton = new wxButton(panel, wxID_ANY, "Rozwiąż", wxDefaultPosition, wxSize(200, 50));
    solveButton->Bind(wxEVT_BUTTON, &JacobiFrame::OnSolve, this);
    
    wxBoxSizer* templateSizer = new wxBoxSizer(wxHORIZONTAL);
    templateSizer->Add(templateChoice, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    templateSizer->Add(templateButton, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    
    rightSizer->Add(templateSizer, 0, wxALL, 5);
    rightSizer->Add(modeRadio, 0, wxALL, 5);
    rightSizer->Add(solveButton, 0, wxALL, 5);
    
    topSizer->Add(inputSizer, 1, wxEXPAND);
    topSizer->Add(rightSizer, 0, wxEXPAND);
    
    mainSizer->Add(topSizer, 0, wxEXPAND);
    mainSizer->Add(new wxStaticText(panel, wxID_ANY, "Wyniki:"), 0, wxALL, 5);
    
    outputMemo = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(-1, 300), wxTE_MULTILINE | wxTE_READONLY);
    
    // Ustaw czcionkę o stałej szerokości dla lepszego wyrównania
    wxFont font = wxSystemSettings::GetFont(wxSYS_ANSI_FIXED_FONT);
    font.SetPointSize(16);
    outputMemo->SetFont(font);
    
    mainSizer->Add(outputMemo, 1, wxEXPAND | wxALL, 5);
    
    panel->SetSizer(mainSizer);
}

string JacobiFrame::Trim(const string& str) {
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first) return str;
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

bool JacobiFrame::ParseInputText(ParsedData& data) {
    string text = inputMemo->GetValue().ToStdString();
    
    // Read base parameters from dedicated inputs
    try { data.n = stoi(nInput->GetValue().ToStdString()); } catch(...) {}
    try { data.mit = stoi(mitInput->GetValue().ToStdString()); } catch(...) {}
    try { data.eps = stold(epsInput->GetValue().ToStdString()); } catch(...) {}
    
    // Overwrite with parameters from text if they exist
    regex r_n("n\\s*=\\s*(\\d+)");
    regex r_mit("mit\\s*=\\s*(\\d+)");
    regex r_eps("eps\\s*=\\s*([\\d\\.eE\\+\\-]+)");
    
    smatch match;
    if (regex_search(text, match, r_n)) {
        data.n = stoi(match[1].str());
        nInput->SetValue(wxString::Format("%d", data.n));
    }
    if (regex_search(text, match, r_mit)) {
        data.mit = stoi(match[1].str());
        mitInput->SetValue(wxString::Format("%d", data.mit));
    }
    if (regex_search(text, match, r_eps)) {
        data.eps = stold(match[1].str());
        epsInput->SetValue(match[1].str());
    }
    
    if (data.n <= 0) {
        wxMessageBox("Wymiar n musi być > 0.");
        return false;
    }
    
    regex r_a("a\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]\\s*=\\s*(\\[[^\\]]+\\]|[^,;\\n\\s]+)");
    regex r_b("b\\[\\s*(\\d+)\\s*\\]\\s*=\\s*(\\[[^\\]]+\\]|[^,;\\n\\s]+)");
    regex r_x("x\\[\\s*(\\d+)\\s*\\]\\s*=\\s*(\\[[^\\]]+\\]|[^,;\\n\\s]+)");
    
    auto a_begin = sregex_iterator(text.begin(), text.end(), r_a);
    auto a_end = sregex_iterator();
    for (sregex_iterator i = a_begin; i != a_end; ++i) {
        smatch m = *i;
        data.a_indices.push_back({stoi(m[1].str()), stoi(m[2].str())});
        data.a_values.push_back(Trim(m[3].str()));
    }
    
    auto b_begin = sregex_iterator(text.begin(), text.end(), r_b);
    for (sregex_iterator i = b_begin; i != a_end; ++i) {
        smatch m = *i;
        data.b_indices.push_back(stoi(m[1].str()));
        data.b_values.push_back(Trim(m[2].str()));
    }
    
    auto x_begin = sregex_iterator(text.begin(), text.end(), r_x);
    for (sregex_iterator i = x_begin; i != a_end; ++i) {
        smatch m = *i;
        data.x_indices.push_back(stoi(m[1].str()));
        data.x_values.push_back(Trim(m[2].str()));
    }
    
    return true;
}

void JacobiFrame::OnTemplate(wxCommandEvent& event) {
    int sel = templateChoice->GetSelection();
    
    if (sel == 0) {
        long n;
        if (!nInput->GetValue().ToLong(&n) || n <= 0) {
            wxMessageBox("Podaj poprawny wymiar n w przygotowanym polu.", "Błąd", wxICON_ERROR);
            return;
        }
        
        ostringstream oss;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                oss << "a[" << i << ", " << j << "] = 0";
                if (j < n) oss << ", ";
            }
            oss << "\n";
        }
        oss << "\n";
        for (int i = 1; i <= n; ++i) {
            oss << "b[" << i << "] = 0";
            if (i < n) oss << ", ";
        }
        oss << "\n\n";
        for (int i = 1; i <= n; ++i) {
            oss << "x[" << i << "] = 0";
            if (i < n) oss << ", ";
        }
        oss << "\n";
        
        inputMemo->SetValue(oss.str());
    } else if (sel == 1) {
        nInput->SetValue("4");
        mitInput->SetValue("100");
        epsInput->SetValue("1e-14");
        string tpl = 
            "a[1, 1] = 0, a[1, 2] = 0, a[1, 3] = 1, a[1, 4] = 2\n"
            "a[2, 1] = 2, a[2, 2] = 1, a[2, 3] = 0, a[2, 4] = 2\n"
            "a[3, 1] = 7, a[3, 2] = 3, a[3, 3] = 0, a[3, 4] = 1\n"
            "a[4, 1] = 0, a[4, 2] = 5, a[4, 3] = 0, a[4, 4] = 0\n\n"
            "b[1] = 1, b[2] = 1, b[3] = 1, b[4] = 1\n\n"
            "x[1] = 0, x[2] = 0, x[3] = 0, x[4] = 0\n";
        inputMemo->SetValue(tpl);
    } else if (sel == 2) {
        nInput->SetValue("4");
        mitInput->SetValue("10");
        epsInput->SetValue("1e-14");
        string tpl = 
            "a[1, 1] = -12.235, a[1, 2] = 1.229, a[1, 3] = 0.5597, a[1, 4] = 0\n"
            "a[2, 1] = 1.229, a[2, 2] = -6.78, a[2, 3] = 0.765, a[2, 4] = 0\n"
            "a[3, 1] = 0.5597, a[3, 2] = 0.765, a[3, 3] = 91.0096, a[3, 4] = 2\n"
            "a[4, 1] = 0, a[4, 2] = 0, a[4, 3] = -2, a[4, 4] = 5.5\n\n"
            "b[1] = 0.956, b[2] = 51.5603, b[3] = 2, b[4] = 5.8\n\n"
            "x[1] = 2, x[2] = 0.75, x[3] = -1, x[4] = 0.9\n";
        inputMemo->SetValue(tpl);
    } else if (sel == 3) {
        nInput->SetValue("4");
        mitInput->SetValue("100");
        epsInput->SetValue("1e-14");
        string tpl = 
            "a[1, 1] = -12.235, a[1, 2] = 1.229, a[1, 3] = 0.5597, a[1, 4] = 0\n"
            "a[2, 1] = 1.229, a[2, 2] = -6.78, a[2, 3] = 0.765, a[2, 4] = 0\n"
            "a[3, 1] = 0.5597, a[3, 2] = 0.765, a[3, 3] = 91.0096, a[3, 4] = 2\n"
            "a[4, 1] = 0, a[4, 2] = 0, a[4, 3] = -2, a[4, 4] = 5.5\n\n"
            "b[1] = 0.956, b[2] = 51.5603, b[3] = 2, b[4] = 5.8\n\n"
            "x[1] = 2, x[2] = 0.75, x[3] = -1, x[4] = 0.9\n";
        inputMemo->SetValue(tpl);
    }
}

void JacobiFrame::OnSolve(wxCommandEvent& event) {
    outputMemo->Clear();
    
    ParsedData data;
    if (!ParseInputText(data)) return;
    
    int n = data.n;
    int mode = modeRadio->GetSelection();
    
    MatrixNormal ANorm(n, vector<long double>(n, 0.0));
    VectorNormal BNorm(n, 0.0);
    VectorNormal XNorm(n, 0.0);
    
    MatrixInterval AInt(n, vector<Interval<long double>>(n, Interval<long double>(0,0)));
    VectorInterval BInt(n, Interval<long double>(0,0));
    VectorInterval XInt(n, Interval<long double>(0,0));
    
    JacobiResult res;
    try {
        if (mode == 0) {
            for (size_t k = 0; k < data.a_indices.size(); ++k) {
                int i = data.a_indices[k].first - 1;
                int j = data.a_indices[k].second - 1;
                if (i >= 0 && i < n && j >= 0 && j < n) {
                    ANorm[i][j] = stold(data.a_values[k]);
                }
            }
            for (size_t k = 0; k < data.b_indices.size(); ++k) {
                int i = data.b_indices[k] - 1;
                if (i >= 0 && i < n) BNorm[i] = stold(data.b_values[k]);
            }
            for (size_t k = 0; k < data.x_indices.size(); ++k) {
                int i = data.x_indices[k] - 1;
                if (i >= 0 && i < n) XNorm[i] = stold(data.x_values[k]);
            }
            
            res = SolveJacobiNormal(ANorm, BNorm, XNorm, data.mit, data.eps);
            
        } else if (mode == 1 || mode == 2) {
            auto parseIntervalString = [this](const string& valStr) -> Interval<long double> {
                regex r_iv("\\[\\s*([^,]+)\\s*,\\s*([^\\]]+)\\s*\\]");
                smatch m;
                if (regex_search(valStr, m, r_iv)) {
                    Interval<long double> iv;
                    iv.a = LeftRead<long double>(Trim(m[1].str()));
                    iv.b = RightRead<long double>(Trim(m[2].str()));
                    return iv;
                }
                return IntRead<long double>(valStr);
            };
            
            for (size_t k = 0; k < data.a_indices.size(); ++k) {
                int i = data.a_indices[k].first - 1;
                int j = data.a_indices[k].second - 1;
                if (i >= 0 && i < n && j >= 0 && j < n) {
                    AInt[i][j] = parseIntervalString(data.a_values[k]);
                }
            }
            for (size_t k = 0; k < data.b_indices.size(); ++k) {
                int i = data.b_indices[k] - 1;
                if (i >= 0 && i < n) BInt[i] = parseIntervalString(data.b_values[k]);
            }
            for (size_t k = 0; k < data.x_indices.size(); ++k) {
                int i = data.x_indices[k] - 1;
                if (i >= 0 && i < n) XInt[i] = parseIntervalString(data.x_values[k]);
            }
            
            res = SolveJacobiInterval(AInt, BInt, XInt, data.mit, data.eps);
        }
        
        // Final result format
        if (res.st == 1 || res.st == 2) {
            outputMemo->AppendText(res.message + "\n");
        } else {
            string finalOut = "";
            for (size_t i = 0; i < res.finalX.size(); ++i) {
                finalOut += "x[" + to_string(i+1) + "] = " + res.finalX[i] + "\n";
            }
            finalOut += "it = " + to_string(res.it) + "\nst = " + to_string(res.st) + "\n";
            outputMemo->AppendText(finalOut);
        }
        
    } catch (const exception& e) {
        wxMessageBox(wxString("Wystąpił błąd parsowania/obliczeń: ") + e.what(), "Błąd", wxICON_ERROR);
    }
}
