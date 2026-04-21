# Tratamiento Digital de Imágenes (TDI) 📸

Este repositorio contiene una suite de herramientas de procesamiento visual desarrolladas en **C++**. El enfoque principal es el uso de algoritmos de **bajo nivel de abstracción**, permitiendo comprender la lógica interna del tratamiento de datos antes de pasar a frameworks de alto nivel.

## 🚀 Funcionalidades Principales

### Análisis Morfológico y Conteo
El núcleo del proyecto es un sistema de **segmentación y etiquetado de componentes conexas**, optimizado para el recuento de objetos basado en:
- Umbralización manual y adaptativa.
- Resolución de tablas de equivalencias en dos pasadas.
- Filtrado por área mínima para eliminación de ruido residual.

### Procesamiento de Señal Visual
- **Suavizado:** Filtros de media para reducción de ruido gaussiano.
- **Orden Estadístico:** Filtro de mediana para eliminación de ruido sal y pimienta.
- **Bordes:** Extracción de gradientes mediante el operador de Sobel.
- **Morfología:** Operaciones de erosión y dilatación (máximos).

## 🛠️ Aspectos Técnicos

- **Lenguaje:** C++ (Uso de `std::filesystem` para gestión de directorios).
- **Abstracción:** Implementación manual de kernels y operaciones matriciales.
- **Organización:** Separación automática de archivos en `Resultados`, `Pasos_intermedios` y `Otros_algoritmos`.

## 📁 Instalación y Uso

1. Clone el repositorio.
2. Ejecute **TDI.cpp** usando Visual Studio 2026.
