# Proyecto de Tratamiento Digital de Imágenes (TDI)

## Universidad de Almería - Grado en Ingeniería Informática

Este repositorio contiene la implementación de un sistema de procesamiento de imágenes orientado al conteo y análisis de objetos (cartas), desarrollado como proyecto práctico para la asignatura de **Tratamiento Digital de Imágenes**.

---

## 🚀 Características Principales

- **Umbralizado Inteligente**: Implementación del algoritmo de **Otsu** para binarización automática basada en el análisis del histograma bimodal.
- **Detección de Bordes Avanzada**: Operador de **Sobel Generalizado** implementado mediante la **Pirámide de Pascal**, permitiendo tamaños de máscara configurables ($3 \times 3, 5 \times 5, 7 \times 7$, etc.) para mayor robustez ante el ruido.
- **Etiquetado de Componentes Conexos**: Algoritmo de dos pasos para la identificación, etiquetado y conteo de objetos individuales en imágenes binarizadas.
- **Filtrado por Área**: Discriminación de ruido mediante un umbral de compromiso de **5000 píxeles**, calculado tras un análisis de la invarianza a la escala.
- **Interfaz de Usuario**: Menú interactivo en consola con gestión de errores y acceso directo a archivos mediante la API de Windows.

---

## 🛠️ Estructura del Proyecto

- `/Source`: Código fuente en C++ (`TDI.cpp`, `C_Matrix.h`, `C_Image_BMP.h`).
- `/Run`: Directorio de recursos que contiene las imágenes de prueba (`.bmp`) y donde se generan los resultados del procesamiento.

---

## 📋 Requisitos e Instalación

El proyecto es compatible con entornos Windows (desde Windows 7 hasta Windows 11) usando el IDE **Visual Studio**.

### 1. Clonado del repositorio

```bash
git clone https://github.com/dgp336/Tratamiento-Digital-de-Imagenes.git
```

### 2. Compilación

Se recomienda usar Visual Studio 2019 o superior. Abrir el archivo `TDI.sln` y compilar el proyecto.

---

## 🖥️ Uso del Programa

Al ejecutar la aplicación, se despliega un menú interactivo en consola que permite al usuario gestionar el flujo de trabajo de manera eficiente. El sistema está diseñado para operar de la siguiente forma:

- **Procesamiento de Imágenes (Opciones 1-8)**: 
  Seleccionan automáticamente imágenes de prueba predefinidas (cartas, dados, etc.). El programa ejecuta una cadena de procesamiento secuencial:
  1. **Carga**: Lectura del archivo `.bmp` desde el directorio `/Run`.
  2. **Preprocesamiento**: Conversión a escala de grises y aplicación de un filtro de media para reducción de ruido.
  3. **Segmentación**: Cálculo del umbral óptimo mediante el algoritmo de **Otsu** y binarización.
  4. **Análisis**: Etiquetado de componentes conexos y conteo de objetos que superen el umbral de área definido.
  5. **Salida**: Generación de archivos de imagen resultantes con los objetos identificados.

- **Gestión de Archivos (Opción 9)**: 
  Mediante el uso de la API de Windows (`_WIN32`), esta opción invoca al explorador de archivos para abrir directamente la carpeta `/Run`. Esto facilita al usuario la verificación inmediata de las imágenes de salida sin necesidad de navegar manualmente por el sistema de archivos.

---

## 📚 Fundamentos Técnicos

### 1. Sobel Generalizado mediante la Pirámide de Pascal

Para la detección de bordes, se ha implementado una versión generalizada del operador de Sobel. En lugar de utilizar una máscara fija de $3 \times 3$, el algoritmo genera coeficientes dinámicos basados en la **Pirámide de Pascal**. 
- El vector de suavizado se obtiene de la fila $N-1$ de Pascal, proporcionando una aproximación binomial a un filtro Gaussiano.
- El vector de diferencia se calcula como la derivada discreta de la fila anterior.

Este enfoque permite utilizar máscaras de mayor tamaño (ej. $5 \times 5$ o $7 \times 7$), aumentando significativamente la robustez del sistema frente al ruido de alta frecuencia.

### 2. Algoritmo de Otsu

La binarización no utiliza un umbral fijo, sino que analiza el histograma de la imagen para encontrar el valor que maximiza la varianza entre clases (fondo y objeto). Esto garantiza que el sistema sea capaz de adaptarse a diferentes condiciones de iluminación sin intervención manual.

### 3. Selección del Umbral de Área Crítica

Tras un análisis estadístico sobre la falta de invarianza a la escala del descriptor de área, se ha establecido un valor de **5000 píxeles** como umbral de corte. Este valor se obtuvo mediante la media aritmética de los límites operativos observados:

$$Area_{defecto} = \frac{1749 + 7692}{2} \approx 5000$$

---

## ✒️ Autor

* **David Granados Pérez**
* Grado en Ingeniería Informática, Universidad de Almería (2026).

---

## 📄 Referencias

### Bibliografía

* **Gonzalez, R. C., & Woods, R. E. (2018).** *Digital Image Processing* (4th ed.). Pearson.

### Material Académico

* **Universidad de Almería.** *Diapositivas de la asignatura Tratamiento Digital de Imágenes*. Grado en Ingeniería Informática.

### Recursos Audiovisuales

* **Becker, A. (2017).** *Intro2Robotics: Connected Components in a Binary Image*. [YouTube].
* **Jian Wei Tay. (2020).** *Otsu's Method*. [YouTube].
* **Machine Scribbler. (2020).** *Connected Component Using 2 Pass Algorithm*. [YouTube].
* **The Organic Chemistry Tutor. (2017).** *Standard Deviation Formula and Variance*. [YouTube].
