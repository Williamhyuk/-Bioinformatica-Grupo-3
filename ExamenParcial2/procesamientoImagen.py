import cv2
import os
import numpy as np


#Practica 4
def operador_logaritmico(folder, filename):
    full_filename = os.path.join(folder, filename)
    imagen_resultado = cv2.imread(full_filename)
    imagen_gray = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)

    c = 255 / np.log10(1 + np.max(imagen_gray))
    alto, ancho = imagen_gray.shape

    for x in range(alto):
        for y in range(ancho):
            imagen_resultado[x][y] = c * (np.log10(1 + imagen_gray[x][y]))

    imagen = imagen_resultado
    full_filename_new = os.path.join(folder, 'operloga' + filename)
    cv2.imwrite(full_filename_new, imagen)

    return full_filename_new


def operador_raiz(folder, filename):
    full_filename = os.path.join(folder, filename)
    imagen_resultado = cv2.imread(full_filename)
    imagen_gray = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)

    c = 50 / np.log10(1 + np.max(imagen_gray))
    alto, ancho = imagen_gray.shape

    for i in range(alto):
        for j in range(ancho):
            raiz = c * (np.sqrt(imagen_gray[i][j]))
            if (raiz < 0):  # Si el resultado del pixel es un valor menor que 0
                imagen_resultado[i][j] = 0  # Se le asignara 0
            elif (raiz > 255):  # Si el resultado del pixel es un valor menor que 255
                imagen_resultado[i][j] = 255  # Se le asignara 0
            else:
                imagen_resultado[i][
                    j] = raiz  # Si los valores estan entre el rango de [0,255] se guardan sin mofiicacion

    img_result = imagen_resultado
    full_filename_new = os.path.join(folder, 'operRaiz' + filename)
    cv2.imwrite(full_filename_new, img_result)

    return full_filename_new

#PRACTICA 5
def operador_exponencial(folder, filename, c):
    full_filename = os.path.join(folder, filename)
    imagen_resultado = cv2.imread(full_filename)
    imagen_gray = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)
    b = 1.01

    alto, ancho = imagen_gray.shape

    for x in range(alto):
        for y in range(ancho):
            resultado = c * (np.power(b, imagen_gray[x][y]) - 1)
            if resultado > 255:
                imagen_resultado[x][y] = 255
            else:
                imagen_resultado[x][y] = c * (np.power(b, imagen_gray[x][y]) - 1)

    img_result = imagen_resultado
    full_filename_new = os.path.join(folder, 'operExponen' + filename)
    cv2.imwrite(full_filename_new, img_result)

    return full_filename_new

#pRACTICA 4
def histogram_equalization(folder, filename):
    full_filename = os.path.join(folder, filename)
    imagen_gray = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)

    cantidad_pixeles = imagen_gray.size
    shape = imagen_gray.shape
    height = shape[0]
    width = shape[1]

    L = 256
    S_n = []
    imagen_array1D = imagen_gray.flatten().tolist()

    suma_acumulada = 0

    # Realizamos S_n
    for index in range(L):
        P_n = imagen_array1D.count(index) / cantidad_pixeles
        suma_acumulada = suma_acumulada + P_n
        s_k = int(round(suma_acumulada * (L - 1)))
        S_n.append(s_k)

    # Realizamos el mapeo lineal
    for index in range(cantidad_pixeles):
        imagen_array1D[index] = S_n[imagen_array1D[index]]

    img_result = np.asarray(imagen_array1D)
    img_result = img_result.reshape(height, width)
    full_filename_new = os.path.join(folder, 'HistEqua' + filename)
    cv2.imwrite(full_filename_new, img_result)

    return full_filename_new


def constrast_streching(folder, filename):
    full_filename = os.path.join(folder, filename)  # importante
    imagen_original = cv2.imread(full_filename)  # importante
    imagen_resultado = cv2.imread(full_filename)
    # Detallamos los valores de las variables de Contrast stretching
    a = 0  # límite inferior
    b = 255  # límite superior
    c = np.min(imagen_original)  # El menor valor de los pixeles
    d = np.max(imagen_original)  # El mayor valor de los pixeles

    alto, ancho, canales = imagen_original.shape

    def point_operator(pixel_RGB):
        return (pixel_RGB - c) * ((b - a) / (d - c) + a)

    for x in range(alto):
        for y in range(ancho):
            imagen_resultado[x][y] = point_operator(imagen_original[x][y])

    img_result = imagen_resultado  # importante
    full_filename_new = os.path.join(folder, 'contrastre' + filename)  # importante
    cv2.imwrite(full_filename_new, img_result)  # importante

    return full_filename_new  # importante


def constrast_streching_out(folder, filename):
    full_filename = os.path.join(folder, filename)  # importante
    res = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)  # importante
    imagen_original = cv2.imread(full_filename)
    # Detallamos los valores de las variables de Contrast stretching
    a = 0  # límite inferior
    b = 255  # límite superior
    # c = 69
    # d = 149
    c = np.min(imagen_original)  # El menor valor de los pixeles
    d = np.max(imagen_original)  # El mayor valor de los pixeles

    longi = d - c  # calculamos la longitud del rango
    limitec = (longi * 5) / 100  # calculamos el limite a partir del porcentaje

    newc = c - limitec  # El menor valor  en un limite de 5%

    limited = (longi * 5) / 100  # calculamos el limite a partir del porcentaje
    newd = d + limited

    alto, ancho, canales = imagen_original.shape

    for x in range(alto):
        for y in range(ancho):
            re = (res[x][y] - newc) * ((b - a) / (newd - newc) + a)
            if (re < 0):
                res[x][y] = 0
            elif (re > 255):
                res[x][y] = 255
            else:
                res[x][y] = re

    img_result = res  # importante
    full_filename_new = os.path.join(folder, 'constratreOutlier' + filename)  # importante
    cv2.imwrite(full_filename_new, img_result)  # importante

    return full_filename_new  # importante

#Practica 5
def power_raise(folder, filename, c):
    full_filename = os.path.join(folder, filename)
    imagen_resultado = cv2.imread(full_filename)
    imagen_gray = cv2.imread(full_filename, cv2.IMREAD_GRAYSCALE)

    r = 0.8
    # c = 4
    #c = 2.5
    # #valores de c 1 1.2 0.5 0.8 3 cuando r = 1
    alto, ancho = imagen_gray.shape

    for x in range(alto):
        for y in range(ancho):
            resultado = (c * (np.power(imagen_gray[x][y], r)))
            if resultado > 255:
                imagen_resultado[x][y] = 255
            elif resultado < 0:
                imagen_resultado[x][y] = 0
            else:
                imagen_resultado[x][y] = (c * (np.power(imagen_gray[x][y], r)))

    img_result = imagen_resultado  # importante
    full_filename_new = os.path.join(folder, 'powerRaise' + filename)  # importante
    cv2.imwrite(full_filename_new, img_result)  # importante

    return full_filename_new  # importante

#practica1
def threholding(folder, filename):
    full_filename = os.path.join(folder, filename)
    imagen = cv2.imread(full_filename)
    resultado = cv2.imread(full_filename)

    height, width, chanels = imagen.shape
    limi = 200  # limite inicial
    limf = 230  # limite final

    for i in range(height):
        for j in range(width):
            if (imagen[i][j][0] > limi or imagen[i][j][1] > limf or imagen[i][j][2] < limi):
                resultado[i][j] = 0

    img_result = resultado  # importante
    full_filename_new = os.path.join(folder, 'thres' + filename)  # importante
    cv2.imwrite(full_filename_new, img_result)  # importante

    return full_filename_new  # importante
