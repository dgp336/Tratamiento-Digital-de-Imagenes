#include "C_Image.hpp"
#include "C_Matrix.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <filesystem>
#include <limits>

using namespace std;

// Parámetros por defecto
int areaMinima = 2000;
int kernelMedia = 3;
double dispersionSeguridad = 20.0; //Valor de dispersión a partir del cual la imagen se considera bimodal, para saber si umbralizar o no 

// --- FUNCIONES DE PROCESAMIENTO ---

void NormalizarHistograma(vector<double>& histograma) {
    const int NIVELES = 256;
    double suma = 0.0;
    for (int k = 0; k < NIVELES; k++) {
        suma += histograma[k];
    }
    
    if (suma > 0.0) {
        for (int k = 0; k < NIVELES; k++) {
            histograma[k] /= suma;
        }
    }
}

vector<double> ComputarHistograma(C_Image& img) {
    const int NIVELES = 256;
    vector<double> histograma(NIVELES, 0.0);

    // Contar píxeles por intensidad
    for (C_Image::IndexT fila = img.FirstRow(); fila <= img.LastRow(); fila++) {
        for (C_Image::IndexT columna = img.FirstCol(); columna <= img.LastCol(); columna++) {
            int intensidad = static_cast<int>(lround(img(fila, columna)));
            intensidad = max(0, min(255, intensidad));
            histograma[intensidad] += 1.0;
        }
    }

    return histograma;
}

void Umbralizar(C_Image& img, int umbral) {
    for (C_Image::IndexT fila = img.FirstRow(); fila <= img.LastRow(); fila++) {
        for (C_Image::IndexT columna = img.FirstCol(); columna <= img.LastCol(); columna++) {
            img(fila, columna) = (img(fila, columna) < umbral) ? 0 : 255;
        }
    }
}

int UmbralizarOtsu(C_Image& img) {
    const int NIVELES = 256;
    
    // 1) Histograma normalizado
    vector<double> histograma = ComputarHistograma(img);

    // Normalizar histograma
    NormalizarHistograma(histograma);

    // 2) Sumas acumuladas P(k)
    vector<double> P(NIVELES, 0.0);
    P[0] = histograma[0];
    for (int k = 1; k < NIVELES; k++) {
        P[k] = P[k - 1] + histograma[k];
    }

    // 3) Medias acumuladas m(k)
    vector<double> m(NIVELES, 0.0);
    m[0] = 0.0 * histograma[0];
    for (int k = 1; k < NIVELES; k++) {
        m[k] = m[k - 1] + static_cast<double>(k) * histograma[k];
    }

    // Cálculo acumulado de segundo momento para varianzas
    vector<double> m2(NIVELES, 0.0);
    m2[0] = 0.0;
    for (int k = 1; k < NIVELES; k++) {
        m2[k] = m2[k - 1] + static_cast<double>(k * k) * histograma[k];
    }

    // 4) Media global
    const double mediaGlobal = m[NIVELES - 1];
    const double segundoMomentoGlobal = m2[NIVELES - 1];

    // 5) Vector de varianzas intra-clase
    vector<double> varIntraClase(NIVELES, numeric_limits<double>::infinity());
    int umbralOtsu = 0;
    double minimo = numeric_limits<double>::infinity();
    const double EPS = 1e-12;

    for (int k = 0; k < NIVELES - 1; k++) {
        double w0 = P[k];
        double w1 = 1.0 - w0;

        if (w0 <= EPS || w1 <= EPS) {
            continue;
        }

        double mu0 = m[k] / w0;
        double mu1 = (mediaGlobal - m[k]) / w1;

        double var0 = (m2[k] / w0) - (mu0 * mu0);
        double var1 = ((segundoMomentoGlobal - m2[k]) / w1) - (mu1 * mu1);

        if (var0 < 0.0) var0 = 0.0;
        if (var1 < 0.0) var1 = 0.0;

        varIntraClase[k] = (w0 * var0) + (w1 * var1);

        // 6) Umbral Otsu: k con varianza intra-clase mínima
        if (varIntraClase[k] < minimo) {
            minimo = varIntraClase[k];
            umbralOtsu = k;
        }
    }

    return umbralOtsu;
}

double calcularDesviacionEstandar(C_Image& img) {
    double suma = 0, sumaCuad = 0;
    int n = img.RowN() * img.ColN();
    
    for (C_Image::IndexT r = img.FirstRow(); r <= img.LastRow(); r++) {
        for (C_Image::IndexT c = img.FirstCol(); c <= img.LastCol(); c++) {
            suma += img(r, c);
            sumaCuad += pow(img(r, c), 2);
        }
    }
    double media = suma / n;
    return sqrt((sumaCuad / n) - pow(media, 2));
}

int EtiquetarYFiltrar(C_Image& binaria, int areaMinima) {
    C_Image etiquetas(binaria.FirstRow(), binaria.LastRow(),
        binaria.FirstCol(), binaria.LastCol(), 0);

    const int MAX_LABELS = 5000;
    vector<int> padre(MAX_LABELS);
    for (int i = 0; i < MAX_LABELS; ++i) padre[i] = i;

    int siguienteEtiqueta = 1;

    // --- PRIMERA PASADA (Detección y Unión) ---
    for (int r = binaria.FirstRow(); r <= binaria.LastRow(); r++) {
        for (int c = binaria.FirstCol(); c <= binaria.LastCol(); c++) {
            if (binaria(r, c) == 255) {
                int arriba = (r > binaria.FirstRow()) ? etiquetas(r - 1, c) : 0;
                int izq = (c > binaria.FirstCol()) ? etiquetas(r, c - 1) : 0;

                if (arriba == 0 && izq == 0) {
                    if (siguienteEtiqueta < MAX_LABELS)
                        etiquetas(r, c) = siguienteEtiqueta++;
                }
                else if (arriba != 0 && izq == 0) {
                    etiquetas(r, c) = arriba;
                }
                else if (arriba == 0 && izq != 0) {
                    etiquetas(r, c) = izq;
                }
                else {
                    etiquetas(r, c) = (arriba < izq) ? arriba : izq;
                    if (arriba != izq) {
                        int raizA = arriba; while (padre[raizA] != raizA) raizA = padre[raizA];
                        int raizI = izq;    while (padre[raizI] != raizI) raizI = padre[raizI];
                        if (raizA != raizI) padre[raizI] = raizA;
                    }
                }
            }
        }
    }

    // --- PASO INTERMEDIO: Conteo de Áreas ---
    vector<int> areas(siguienteEtiqueta, 0);
    for (int r = etiquetas.FirstRow(); r <= etiquetas.LastRow(); r++) {
        for (int c = etiquetas.FirstCol(); c <= etiquetas.LastCol(); c++) {
            int e = etiquetas(r, c);
            if (e != 0) {
                int raiz = e;
                while (padre[raiz] != raiz) {
                    padre[raiz] = padre[padre[raiz]];
                    raiz = padre[raiz];
                }
                etiquetas(r, c) = raiz;
                areas[raiz]++;
            }
        }
    }

    // --- SEGUNDA PASADA: Filtrado y Reetiquetado Consecutivo ---
    vector<int> nuevoID(siguienteEtiqueta, 0);
    int contadorReal = 1;

    for (int i = 1; i < siguienteEtiqueta; i++) {
        if (padre[i] == i) {
            if (areas[i] >= areaMinima) {
                nuevoID[i] = contadorReal++;
            }
            else {
                nuevoID[i] = 0;
            }
        }
    }

    for (int r = etiquetas.FirstRow(); r <= etiquetas.LastRow(); r++) {
        for (int c = etiquetas.FirstCol(); c <= etiquetas.LastCol(); c++) {
            int e = etiquetas(r, c);
            if (e != 0) {
                etiquetas(r, c) = nuevoID[e];
            }
        }
    }

    return contadorReal - 1;
}

// --- OTRAS FUNCIONES IMPLEMENTADAS ---

void FiltroMaximos(C_Image& img, int longitudMascara) {
    int offset = longitudMascara / 2;
    C_Image temp(img);

    for (C_Image::IndexT fila = img.FirstRow() + offset; fila <= img.LastRow() - offset; fila++) {
        for (C_Image::IndexT columna = img.FirstCol() + offset; columna <= img.LastCol() - offset; columna++) {
            double maximo = -1.0;
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    double val = temp(fila + i, columna + j);
                    if (val > maximo) {
                        maximo = val;
                    }
                }
            }
            img(fila, columna) = maximo;
        }
    }
}

void FiltroMinimos(C_Image& img, int tamVentana) {
    int offset = tamVentana / 2;
    C_Image temp(img);

    for (C_Image::IndexT fila = img.FirstRow() + offset; fila <= img.LastRow() - offset; fila++) {
        for (C_Image::IndexT columna = img.FirstCol() + offset; columna <= img.LastCol() - offset; columna++) {
            double minimo = 255.0;
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    double valor = temp(fila + i, columna + j);
                    if (valor < minimo) minimo = valor;
                }
            }
            img(fila, columna) = minimo;
        }
    }
}

void FiltroMediana(C_Image& img, int mascara) {
    int offset = mascara / 2;
    int size = mascara * mascara;
    C_Image temp(img);
    vector<double> ventana(size); // Vector preasignado

    for (C_Image::IndexT fila = img.FirstRow() + offset; fila <= img.LastRow() - offset; fila++) {
        for (C_Image::IndexT columna = img.FirstCol() + offset; columna <= img.LastCol() - offset; columna++) {
            int idx = 0;
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    ventana[idx++] = temp(fila + i, columna + j);
                }
            }
            // nth_element es mucho más rápido que sort para encontrar la mediana
            nth_element(ventana.begin(), ventana.begin() + size / 2, ventana.end());
            img(fila, columna) = ventana[size / 2];
        }
    }
}

void DeteccionBordesSobel(C_Image& entrada, C_Image& salida) {
    int sobelX[3][3] = { {-1,0,1}, {-2,0,2}, {-1,0,1} };
    int sobelY[3][3] = { {-1,-2,-1}, {0,0,0}, {1,2,1} };

    salida.Resize(entrada.FirstRow(), entrada.LastRow(), entrada.FirstCol(), entrada.LastCol());

    for (C_Image::IndexT fila = entrada.FirstRow() + 1; fila <= entrada.LastRow() - 1; fila++) {
        for (C_Image::IndexT columna = entrada.FirstCol() + 1; columna <= entrada.LastCol() - 1; columna++) {
            double gx = 0, gy = 0;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    gx += entrada(fila + i, columna + j) * sobelX[i + 1][j + 1];
                    gy += entrada(fila + i, columna + j) * sobelY[i + 1][j + 1];
                }
            }
            salida(fila, columna) = min<double>(255.0, sqrt(gx * gx + gy * gy));
        }
    }
}

// --- FUNCIÓN GENÉRICA DE PROCESAMIENTO PARAMETRIZADA ---

void ProcesarImagen(string nombreArchivo, int areaMinima, int kernelMedia) {
    C_Image img;
    int numCartas;
    string nombreBase = nombreArchivo.substr(0, nombreArchivo.find_last_of("."));

    C_Trace(("Comenzando procedimiento para " + nombreArchivo + "...").c_str());
    img.Read(nombreArchivo.c_str());

    img.Grey();
    C_Trace("Imagen convertida a escala de grises.");
    img.Write(("Pasos intermedios/" + nombreBase + "_01_Gris.bmp").c_str());

    C_Trace(("Aplicando filtro de la media " + to_string(kernelMedia) + "x" + to_string(kernelMedia) + "...").c_str());
    img.MeanFilter(img, kernelMedia);
    img.Write(("Pasos intermedios/" + nombreBase + "_02_Media.bmp").c_str());

    C_Trace("Analizando histograma...");
    double dispersion = calcularDesviacionEstandar(img);

    int umbralFinal;
    if (dispersion < dispersionSeguridad) { // Umbral de seguridad para imágenes "planas"
        C_Trace("Histograma compacto detectado (posible imagen de fondo). Forzando umbral alto.");
        umbralFinal = 250; 
    } else {
        C_Trace("Histograma bimodal detectado. Umbralizando imagen automaticamente con Otsu...");
        umbralFinal = UmbralizarOtsu(img);
        C_Trace(("Umbral Otsu calculado: " + to_string(umbralFinal)).c_str());
    }

    // Ahora aplicas la binarización con el umbral seleccionado
    Umbralizar(img, umbralFinal);
    img.Write(("Resultados/" + nombreBase + "_Procesada_Binaria.bmp").c_str());

    C_Trace(("Inicio del algoritmo de recuento de area mayor que " + to_string(areaMinima) + "...").c_str());
    numCartas = EtiquetarYFiltrar(img, areaMinima);

    C_Trace("\n----------------------------------------\n");
    C_Trace((" RESULTADO PARA " + nombreArchivo).c_str());
    C_Trace((" Cartas detectadas: " + to_string(numCartas)).c_str());
    C_Trace("----------------------------------------\n");

    img.Free();
}

// --- MENÚ DE OTRAS OPERACIONES ---
void MenuOtrasOperaciones() {
    string archivo;
    int opcion;
    int tamanoMascara = 3;

    cout << "\n[OTRAS OPERACIONES DE PROCESAMIENTO]" << endl;

    cout << "\nSeleccione el filtro a aplicar:" << endl;
    cout << "1. Filtro Erosion" << endl;
    cout << "2. Filtro Maximos" << endl;
    cout << "3. Filtro Mediana" << endl;
    cout << "4. Deteccion de Bordes (Sobel)" << endl;
    cout << "0. Volver" << endl;
    cout << "Seleccione: ";
    cin >> opcion;

    if (opcion == 0) return;

    cout << "Introduzca el nombre de la imagen a procesar: ";
    cin >> archivo;

    C_Image img;
    img.Read(archivo.c_str());
    img.Grey(); // Pasamos a gris por defecto para estos filtros
    string nombreBase = archivo.substr(0, archivo.find_last_of("."));

    if (opcion >= 1 && opcion <= 3) {
        cout << "Introduzca el tamano de la mascara (Ej: 3 para 3x3): ";
        cin >> tamanoMascara;
        if (tamanoMascara % 2 == 0) tamanoMascara++; // Forzamos impar
    }

    switch (opcion) {
    case 1:
        FiltroMinimos(img, tamanoMascara);
        img.Write(("Otros algoritmos/" + nombreBase + "_Minimos.bmp").c_str());
        cout << "-> Guardado en 'Otros algoritmos/" << nombreBase << "_Minimos.bmp'" << endl;
        break;
    case 2:
        FiltroMaximos(img, tamanoMascara);
        img.Write(("Otros algoritmos/" + nombreBase + "_Maximos.bmp").c_str());
        cout << "-> Guardado en 'Otros algoritmos/" << nombreBase << "_Maximos.bmp'" << endl;
        break;
    case 3:
        FiltroMediana(img, tamanoMascara);
        img.Write(("Otros algoritmos/" + nombreBase + "_Mediana.bmp").c_str());
        cout << "-> Guardado en 'Otros algoritmos/" << nombreBase << "_Mediana.bmp'" << endl;
        break;
    case 4:
    {
        C_Image salidaSobel;
        DeteccionBordesSobel(img, salidaSobel);
        salidaSobel.Write(("Otros algoritmos/" + nombreBase + "_Sobel.bmp").c_str());
        cout << "-> Guardado en 'Otros algoritmos/" << nombreBase << "_Sobel.bmp'" << endl;
        break;
    }
    default:
        cout << "Opcion no valida." << endl;
        break;
    }
}

// --- FUNCIÓN MAIN ---

int main() {

    // Creación dinámica de carpetas al iniciar el programa
    filesystem::create_directories("Resultados");
    filesystem::create_directories("Pasos intermedios");
    filesystem::create_directories("Otros algoritmos");

    int opcion = 0;

    cout << "========================================" << endl;
    cout << "   SISTEMA DE CONTEO DE CARTAS (TDI)    " << endl;
    cout << "========================================" << endl;

    do {
        cout << "\n----------- PARAMETROS CONFIGURADOS -----------" << endl;
        cout << "\n Area minima de conteo de objeto:" << to_string(areaMinima) << endl;
        cout << "\n Tamano del kernel de la media:" << to_string(kernelMedia) << endl;

        cout << "\n----------- MENU DE OPCIONES -----------" << endl;
        cout << "1. Procesar cartas01.bmp (8 cartas)" << endl;
        cout << "2. Procesar cartas02.bmp (5 cartas)" << endl;
        cout << "3. Procesar cartas03.bmp (6 cartas)" << endl;
        cout << "4. Procesar cartas04.bmp (Imagen de control, 0 cartas)" << endl;
        cout << "5. Procesar todas las imagenes" << endl;
        cout << "6. Configurar parametros" << endl;
        cout << "7. Otras operaciones (Filtros extra)" << endl;
        cout << "0. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
        case 1: ProcesarImagen("cartas01.bmp", areaMinima, kernelMedia); break;
        case 2: ProcesarImagen("cartas02.bmp", areaMinima, kernelMedia); break;
        case 3: ProcesarImagen("cartas03.bmp", areaMinima, kernelMedia); break;
        case 4: ProcesarImagen("cartas04.bmp", areaMinima, kernelMedia); break;
        case 5: 
            for (int i = 1; i <= 4;i++) {
                string carta = "cartas0" + to_string(i) + ".bmp";
                ProcesarImagen(carta, areaMinima, kernelMedia);
            }
            break;
        case 6:
            cout << "\n[CONFIGURACION DE PARAMETROS]" << endl;

            cout << "Introduzca el tamano del kernel para el suavizado de la media (Recomendado > 3): ";
            cin >> kernelMedia;

            cout << "Introduzca el area minima para descontar ruido (Recomendado > 2000): ";
            cin >> areaMinima;

            // Validación básica
            if (kernelMedia % 2 == 0) kernelMedia++; // Forzamos a que sea impar
            if (areaMinima <= 0) areaMinima = 50;

            cout << "\n-> Parametros configurados correctamente.\n";
            break;
        case 7:
            MenuOtrasOperaciones();
            break;
        case 0: cout << "Saliendo del programa..." << endl; break;
        default: cout << "Opcion no valida." << endl; break;
        }

    } while (opcion != 0);

    return 0;
}