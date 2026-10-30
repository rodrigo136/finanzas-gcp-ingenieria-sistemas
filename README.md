# FinanzasGCP — Planificación financiera con Gradientes Conjugados Precondicionados
Estudiante: Marca Vargas Ariel Rodrigo

Esta es una aplicación web académica para el desafio METODO DEL GRADIENTE - AREA FINANCIERA en la materia de Métodos Numéricos orientada a la mención de **Ingeniería de Sistemas**.

## 1. Descripción

El proyecto modela la planificación financiera y distribución ponderada de recursos entre cuatro unidades de una organización tecnológica:

1. Desarrollo de software.
2. Infraestructura tecnológica.
3. Seguridad informática.
4. Soporte y mantenimiento.

El problema se formula como un sistema lineal:

\[
Ax=b
\]

La aplicación permite trabajar con tres escenarios: ideal, bajo estrés financiero y mal condicionado. Compara la solución numérica de referencia obtenida por Cholesky, el método de Gradientes Conjugados (GC) y Gradientes Conjugados Precondicionados (GCP) con precondicionador de Jacobi.

> **Importante:** los datos son sintéticos y se utilizan con fines académicos. No son registros financieros reales ni recomendaciones de gasto.

## 2. Funcionalidades

- Presentación del caso financiero y justificación de variables.
- Definición de cuatro variables financieras y explicación de los coeficientes.
- Carga de escenarios predefinidos.
- Edición interactiva de la matriz simétrica `A` y del vector `b`.
- Configuración de vector inicial, tolerancia, precondicionador y máximo de iteraciones.
- Validación de simetría y definición positiva mediante factorización de Cholesky.
- Solución directa de referencia mediante Cholesky.
- Implementación de GC simple y GCP con Jacobi.
- Cálculo de residuo relativo, número de condición y condición precondicionada.
- Tabla de soluciones y errores absoluto, relativo y porcentual.
- Gráfico de convergencia de residuos.
- Módulo QA con cuatro pruebas funcionales.
- Diseño adaptable a escritorio, tableta y móvil.

## 3. Tecnologías

- HTML5: estructura y contenido.
- CSS3: diseño adaptable.
- JavaScript ES6: lógica de interfaz y algoritmos numéricos.
- Canvas 2D: gráfico de convergencia.

No se requiere instalar dependencias ni utilizar un servidor backend.

## 4. Estructura del proyecto

```text
finanzas-gcp-ingenieria-sistemas/
├── index.html
├── README.md
└── src/
    ├── app.js
    ├── math.js
    └── styles.css
```

### Responsabilidad de cada archivo

| Archivo | Responsabilidad |
|---|---|
| `index.html` | Estructura de la página, secciones, tablas, controles y accesibilidad básica. |
| `src/styles.css` | Identidad visual, distribución de componentes y adaptación responsive. |
| `src/math.js` | Operaciones matriciales, Cholesky, autovalores, condición, GC y GCP. |
| `src/app.js` | Escenarios, edición de datos, validaciones, cálculos, resultados, gráficos y pruebas QA. |
| `README.md` | Documentación, instalación, ejecución y explicación del proyecto. |

## 5. Modelo y variables

Las incógnitas representan asignaciones financieras ponderadas, expresadas como miles de bolivianos:

- `x1`: desarrollo de software.
- `x2`: infraestructura tecnológica.
- `x3`: seguridad informática.
- `x4`: soporte y mantenimiento.

Los coeficientes diagonales `aᵢᵢ` ponderan la necesidad y capacidad de absorción de recursos de cada unidad. Los coeficientes fuera de la diagonal `aᵢⱼ` representan interdependencias y costos compartidos. El vector `b` representa metas financieras ponderadas.

Los parámetros no son montos contables directos: constituyen un modelo sintético para demostrar la resolución numérica. En una implementación institucional real, deberían calibrarse con datos históricos, costos, prioridades y reglas de negocio.

### Escenarios

- **Ideal:** matriz simétrica equilibrada y vector de metas estable.
- **Bajo estrés:** se conserva la matriz y se incrementa el vector `b` para representar inflación, incremento de costos tecnológicos y presión operativa.
- **Mal condicionado:** se emplea una matriz `A = D C D`, donde `D` tiene escalas diagonales muy distintas. Esto crea mal escalamiento diagonal y justifica el precondicionador de Jacobi.

Las matrices de los escenarios son simétricas y definidas positivas. La aplicación vuelve a comprobar estas propiedades al ejecutar los cálculos; la edición manual puede producir matrices inválidas, que son rechazadas con un mensaje.

## 6. Métodos numéricos implementados

### Solución de referencia: Cholesky

Para una matriz simétrica definida positiva se factoriza:

\[
A=LL^T
\]

Luego se resuelven los sistemas triangulares:

\[
Ly=b,\qquad L^Tx=y
\]

### Gradientes Conjugados (GC)

Se resuelve el sistema sin precondicionamiento. Se utiliza el residuo inicial `r0 = b - A x0` y se actualiza la solución hasta satisfacer:

\[
\frac{\|r_k\|_2}{\|b\|_2}<\text{tolerancia}
\]

### Gradientes Conjugados Precondicionados (GCP)

El precondicionador de Jacobi utiliza:

\[
M=\operatorname{diag}(A)
\]

En cada iteración se resuelve implícitamente `M z = r` mediante división componente a componente, sin construir la inversa de `M`.

### Número de condición

Se calcula el número de condición espectral para matrices simétricas definidas positivas:

\[
\kappa_2(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}
\]

Para Jacobi, se calcula mediante la matriz simétrica equivalente:

\[
M^{-1/2} A M^{-1/2}
\]

que tiene los mismos autovalores que `M⁻¹A`.

### Errores

Para cada variable se presentan:

- Error absoluto: `Ea = |xref - xGCP|`.
- Error relativo: `Er = Ea / |xref|`.
- Error porcentual: `Ep = Er × 100`.
- Residuo relativo global: `||b - A xGCP||₂ / ||b||₂`.

Si la solución de referencia de una componente es cero, el error relativo de esa componente se informa como no definido.

## 7. Requisitos

- Un navegador moderno (Chrome, Firefox, Edge o Safari).
- Git, si se desea clonar o publicar el repositorio.
- No se requiere Node.js, npm, base de datos ni servidor backend.

## 8. Instalación y ejecución local

### Opción A: descarga del repositorio

```bash
git clone https://github.com/TU-USUARIO/finanzas-gcp-ingenieria-sistemas.git
cd finanzas-gcp-ingenieria-sistemas
```

Abre `index.html` directamente en el navegador.

### Opción B: servidor local opcional

Si tienes Python instalado, desde la carpeta del proyecto ejecuta:

```bash
python -m http.server 8000
```

Abre en el navegador:

```text
http://localhost:8000
```

## 9. Publicación en GitHub Pages

1. Crea un repositorio público en GitHub, por ejemplo `finanzas-gcp-ingenieria-sistemas`.
2. Copia todos los archivos y carpetas del proyecto al repositorio.
3. Realiza el primer commit y push:

```bash
git init
git add .
git commit -m "feat: implementar simulador financiero GCP"
git branch -M main
git remote add origin https://github.com/TU-USUARIO/finanzas-gcp-ingenieria-sistemas.git
git push -u origin main
```

4. En GitHub entra a **Settings → Pages**.
5. En **Build and deployment**, selecciona **Deploy from a branch**.
6. Selecciona la rama `main` y la carpeta `/ (root)`.
7. Guarda los cambios y espera a que GitHub Pages publique el sitio.

## 10. Uso de la aplicación

1. Entra a la sección **Simulador numérico**.
2. Selecciona Ideal, Bajo estrés o Mal condicionado.
3. Revisa la matriz `A` y el vector `b`.
4. Modifica valores si deseas experimentar. La matriz debe continuar siendo simétrica y definida positiva para ejecutar GC/GCP.
5. Configura el vector inicial, tolerancia, precondicionador y máximo de iteraciones.
6. Presiona **Ejecutar cálculo**.
7. Consulta el diagnóstico, gráfico de convergencia, tabla comparativa, errores y explicación financiera.
8. En **Módulo QA**, presiona **Ejecutar pruebas QA** para registrar los resultados de los cuatro casos de prueba.

## 11. Pruebas QA incluidas

| Caso | Entrada | Resultado esperado |
|---|---|---|
| Carga de escenarios | Ideal, estrés y mal condicionado | Cada escenario debe ser válido y resoluble. |
| Matriz no válida | Matriz no simétrica | La aplicación debe detectar y notificar el problema. |
| Edición dinámica | Alterar un coeficiente de `A` | La solución recalculada debe reflejar el cambio. |
| Precondicionamiento | Escenario mal condicionado con Jacobi | GCP debe alcanzar la tolerancia en menor o igual número de iteraciones que GC para el caso de prueba configurado. |

Las pruebas se ejecutan en el navegador y muestran el resultado obtenido en la tabla QA. No sustituyen pruebas unitarias automatizadas en un pipeline CI.

## 12. Limitaciones y consideraciones

- Los datos son hipotéticos y no deben utilizarse para decisiones presupuestarias reales.
- La dimensión del sistema está fijada en cuatro variables para ajustarse al caso de estudio.
- El precondicionador disponible es Jacobi; la opción “Sin precondicionador” permite contrastar con GC.
- El cálculo del número de condición utiliza autovalores de una matriz simétrica mediante rotaciones de Jacobi.
- El método directo Cholesky y los métodos iterativos se implementan en JavaScript para ejecución local en el navegador.
- La solución matemática no incluye restricciones adicionales como presupuesto total, límites mínimos/máximos o variables enteras. Esas restricciones requerirían un modelo de optimización distinto.

## 13. Bibliografía

- Suñagua S., P. (2020). *Método de Gradientes Conjugados Precondicionado*. Revista Boliviana de Matemática – UMSA, 04, 2–7.
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.). SIAM.

---
