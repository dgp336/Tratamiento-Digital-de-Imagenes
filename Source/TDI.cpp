#include "C_Image.hpp"
#include "C_Matrix.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <filesystem>

using namespace std;

// --- FUNCIONES DE PROCESAMIENTO ---

void UmbralizarFondo(C_Image& img, int umbral) {
    for (C_Image::IndexT fila = img.FirstRow(); fila <= img.LastRow(); fila++) {
        for (C_Image::IndexT columna = img.FirstCol(); columna <= img.LastCol(); columna++) {
            img(fila, columna) = (img(fila, columna) < umbral) ? 0 : 255;
        }
    }
}

int etiquetarYFiltrar(C_Image& binaria, int areaMinima) {
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

void filtroMaximos(C_Image& img, int longitudMascara) {
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

void FiltroErosion(C_Image& img, int tamVentana) {
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
            // Optimización: nth_element es mucho más rápido que sort para encontrar la mediana
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

void procesarImagen(string nombreArchivo, int areaMinima, int umbral, int kernelMedia) {
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

    C_Trace(("Umbralizando imagen con valor " + to_string(umbral) + "...").c_str());
    UmbralizarFondo(img, umbral);
    img.Write(("Resultados/" + nombreBase + "_Procesada_Binaria.bmp").c_str());

    C_Trace(("Inicio del algoritmo de recuento de area mayor que" + to_string(areaMinima) + "...").c_str());
    numCartas = etiquetarYFiltrar(img, areaMinima);

    C_Trace("\n----------------------------------------\n");
    C_Trace((" RESULTADO PARA " + nombreArchivo).c_str());
    C_Trace((" Cartas detectadas: " + to_string(numCartas)).c_str());
    C_Trace("----------------------------------------\n");

    img.Free();
}

// --- MENÚ DE OTRAS OPERACIONES ---
void menuOtrasOperaciones() {
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
        FiltroErosion(img, tamanoMascara);
        img.Write(("Otros algoritmos/" + nombreBase + "_Erosion.bmp").c_str());
        cout << "-> Guardado en 'Otros algoritmos/" << nombreBase << "_Erosion.bmp'" << endl;
        break;
    case 2:
        filtroMaximos(img, tamanoMascara);
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

    // Parámetros por defecto
    int areaMinima = 2000;
    int umbral = 90;
    int kernelMedia = 3;

    cout << "========================================" << endl;
    cout << "   SISTEMA DE CONTEO DE CARTAS (TDI)    " << endl;
    cout << "========================================" << endl;

    do {
        cout << "\n----------- PARAMETROS CONFIGURADOS -----------" << endl;
        cout << "\n Umbral de binarizacion:" << to_string(umbral) << endl;
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
        case 1: procesarImagen("cartas01.bmp", areaMinima, umbral, kernelMedia); break;
        case 2: procesarImagen("cartas02.bmp", areaMinima, umbral, kernelMedia); break;
        case 3: procesarImagen("cartas03.bmp", areaMinima, umbral, kernelMedia); break;
        case 4: procesarImagen("cartas04.bmp", areaMinima, umbral, kernelMedia); break;
        case 5: 
            for (int i = 1; i <= 4;i++) {
                string carta = "cartas0" + to_string(i) + ".bmp";
                procesarImagen(carta, areaMinima, umbral, kernelMedia);
            }
            break;
        case 6:
            cout << "\n[CONFIGURACION DE PARAMETROS]" << endl;

            cout << "Introduzca el umbral de binarizacion (Ej: 90): ";
            cin >> umbral;

            cout << "Introduzca el kernel para el suavizado de la media (Recomendado > 5): ";
            cin >> kernelMedia;

            cout << "Introduzca el area minima para descontar ruido (Recomendado > 2000): ";
            cin >> areaMinima;

            // Validación básica
            if (kernelMedia % 2 == 0) kernelMedia++;       // Forzamos a que sea impar
            if (areaMinima <= 0) areaMinima = 50;

            cout << "\n-> Parametros configurados correctamente.\n";
            break;
        case 7:
            menuOtrasOperaciones();
            break;
        case 0: cout << "Saliendo del programa..." << endl; break;
        default: cout << "Opcion no valida." << endl; break;
        }

    } while (opcion != 0);

    return 0;
}