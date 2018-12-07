import numpy as np
import cv2 as cv
from sympy import *
from numpy.linalg import inv


#coordinates of principal point im milimeters
x0 = -0.223
y0 = 0.051
f = 50
var_p = 0.0046
it = 1

#Radial Distortion coeficients

K1 = -1.5604645*10**-5
K2 = 4.0643961*10**-9
K3 = 0

#Decentering Distortion coeficients P1 and P2

P1 = 0
P2 = 0

# Functional Model: Collinearity equations
# Rotation order: R(z) R(y) R(x)
	
X, Y, Z, X0, Y0, Z0, om, fi, kapa, c  = symbols('X Y Z X0 Y0 Z0 om fi kapa c')

m11 = cos(fi) * cos(kapa)
m12 = cos(om) *sin(kapa) + sin(om)*sin(fi)*cos(kapa)
m13 = sin(om) *sin(kapa) - cos(om)*sin(fi)*cos(kapa)

m21 = -cos(fi)* sin(kapa)
m22 = cos(om)*cos(kapa) - sin(om)*sin(fi)*sin(kapa)
m23 = sin(om)*cos(kapa) + cos(om)*sin(fi)*sin(kapa)

m31 = sin(om)
m32 = -sin(om)*cos(fi)
m33 = cos(om)*cos(fi)


x = - c * (m11*(X-X0) + m12*(Y-Y0) + m13*(Z-Z0)) / (m31*(X-X0) + m32*(Y-Y0) + m33*(Z-Z0))
y = - c * (m21*(X-X0) + m22*(Y-Y0) + m23*(Z-Z0)) / (m31*(X-X0) + m32*(Y-Y0) + m33*(Z-Z0))

#---------------------------------- ---------------------------------------------------------------
#----------------------------------Matriz de Pesos ------------------------------------------------
#--------------------------------------------------------------------------------------------------

lin_aux = []
matriz_P = []

#Create weight matrix
for i in range(32):
	for j in range(32):
		if i==j:
			lin_aux.append(var_p/2)
		else:
			lin_aux.append(0)
	matriz_P.append(lin_aux)
	lin_aux = []


def calculateRadialDist(x_obs, y_obs):

	#calculating r
	r = np.sqrt(x_obs**2 + y_obs**2)

	#calculating dx and dy
	dx = (K1*r**2 + K2*r**4 + K3*r**6)*x_obs
	dy = (K1*r**2 + K2*r**4 + K3*r**6)*y_obs

	x_corr = x_obs - dx
	y_corr = y_obs - dy

	corr_radialDist = [x_corr, y_corr]

	return corr_radialDist

def calculateDecentringDist(x_obs, y_obs):

	#calculating r
	r = np.sqrt(x_obs**2 + y_obs**2)

	#calculate dx and dy
	dx = P1*(r**2 + 2*x_obs**2) + 2*P2*x_obs*y_obs
	dy = 2*P1*x_obs*y_obs + P2*(r**2 + 2*y_obs**2)

	decentringDist = [dx, dy]

	return decentringDist



def photogrametricCorrection(radialList, decentringList, x_obs, y_obs): #sei lá por que eu fiz isso aqui

	x_corrected = x_obs + radialList[0] + decentringList[0]
	y_corrected = y_obs + radialList[1] + decentringList[1]

	correctedPhotocoordinates = [x_corrected, y_corrected]

	return correctedPhotocoordinates

def convertCoordinates(C, L):

	img = cv.imread('img1.tif'); #You must give the path of your image

	rows, columms, bands = img.shape #The method shape return the dimensions of the image, but here the bands are not interesting

	#print("The image size: ")
	#print(img.shape)

	pixel_size_x = 0.0046
	pixel_size_y = 0.0046

	x_mm = 0
	y_mm = 0

	x_mm = pixel_size_x * (C - ((columms-1)/2)) #calculate the x coordinate
	y_mm = -pixel_size_y * (L - ((rows-1)/2)) #calculate the y coordinate

	#print("The transformed coordinates: ")
	#print(x_mm, y_mm)

	transformed = [x_mm,y_mm]

	return transformed

def readArqs():

	'''Esta função lê as observações realizadas nas duas fotografias, tanto para os pontos tie quanto para os GCP. Veja a planilha 
	   para saber como montar o arquivo foto_observ.txt
	 '''
	 '''
	 	This function reads the observations realised on the two photos. See the file foto_observ.txt to see how to organize the file.
		 - The first 4 collumns are the C,L of the Tie points in the digital image coordinate system in the photo 1 and photo 2 respectively
		 - The 5th, 6th and 7th collumns are the initial aproximations for X, Y, Z of the tie points in a geodetic reference system
		 - The next 6 collumns is the initial aproximations for the exterior orientations parameters X0, Y0, Z0 of the two photos respectively,
		   i.e., 8th, 9th and 10th collumns is the X0, Y0, Z0 of the first photo, and 11th, 12th and 13th collumns is X0, Y0 and Z0 of the second photo.
		 - The next 6 collumns is the initial aproximations for the attitude angles (omega, phi, kappa) of the first and second photo respectively (similar way of X0,Y0 and Z0)

	 '''

	arq = open('foto_observ.txt', 'r')

	dados1 = arq.readlines() #Armazena o arquivo como um todo

	if dados1:
	   print("Leitura realizada com sucesso!")
	else:
	   print("Leitura falhou!")

	lista_pontos = [] #Lista que armazena cada linha do arquivo como float
	lista_linhas = []

	for i in dados1:
	    dados_separados = i.split(" ") #captura cada numero separado por espaço
	    for j in dados_separados:
	        lista_linhas.append(float(j)) #converte cada número de caractere para float
	        #print(j)
	    lista_pontos.append(lista_linhas) # esta é a lista que interessa
	    lista_linhas = []

	return lista_pontos

#Alocando espaço para a matriz A


#------------------------------------------------------------------------------------------------------------------
#-------------------------------------------Mouting A matrix: Photo 1 -------------------------------------------
#------------------------------------------------------------------------------------------------------------------

def calcula_matA_foto1(parametros):

	a1 = diff(x,X0)
	a2 = diff(x,Y0)
	a3 = diff(x,Z0)
	a4 = diff(x,om)
	a5 = diff(x,fi)
	a6 = diff(x,kapa)
	a7 = diff(y,X0)
	a8 = diff(y,Y0)
	a9 = diff(y,Z0)
	a10 = diff(y,om)
	a11 = diff(y,fi)
	a12 = diff(y,kapa)
	a13 = diff(x, X) #tie
	a14 = diff(x, Y) #tie
	a15 = diff(x, Z) #tie
	a16 = diff(y, X) #tie
	a17 = diff(y, Y) #tie
	a18 = diff(y, Z) #tie

	coef = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18]

	cont = 0 #contador para as linhas do arquivo lido
	cont2 = 0 #Contador para controlar o calculo da lista coef
	val_a = []

	#Calculating the differential coeficients in the parameters

	for i in parametros:
		for j in coef:
			if (cont <= 3):
				num_a = j.evalf(subs={c: 50, X0: i[7], Y0: i[8], Z0: i[9], X: i[4], Y: i[5], Z: i[6], om: i[13], fi: i[14], kapa: i[15]})
				val_a.append(num_a)
				#if cont ==0:
					#print(num_a)
			else:
				if cont2 > 11:
					break
				else:
					num_a = j.evalf(subs={c: 50, X0: i[7], Y0: i[8], Z0: i[9], X: i[4], Y: i[5], Z: i[6], om: i[13], fi: i[14], kapa: i[15]})
					val_a.append(num_a)
					cont2 = cont2 + 1
		cont2 = 0
		cont = cont + 1
	
	
	#print('tamanho coefs: ', len(val_a))
	tie_list = val_a[0:72]
	GCP_list = val_a[72:120]

	#------------------------------------------------------------------------------------------------------------------
	#-------------------------------------------Mouting A matrix: Tie points part ------------------------------------
	#------------------------------------------------------------------------------------------------------------------

	tie_matrix = np.array(tie_list).reshape(12,6)
	tie_matrix_part_list = []
	A_matrix = []
	p1 = np.zeros((2,3), dtype = 'float')
	cont = 0
	cont2 = 0 #controlador da parte dos tie

	tie_matrix_final = []

	for i in tie_matrix:
		if cont != 2:
			tie_matrix_part_list.append(i)
			cont = cont + 1
		else:
			if cont2 == 0:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((tie_matrix_part, p1, p1), axis = 1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p2, p1, p1,p1), axis =1) #funciona adequadamente
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 == 1:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((tie_matrix_part, p1, p1), axis =1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p2, p1, p1), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 ==2:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((tie_matrix_part, p1, p1), axis =1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p1, p2, p1), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 == 3:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((tie_matrix_part, p1, p1), axis =1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p1, p1, p2), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
		
		

	#print(p1)
	#print("--------------------------------------------")
	#print(tie_matrix_part)
	#print("--------------------------------------------")
	#print(p2)
	#print("--------------------------------------------")
	#print(part)
	#print("--------------------------------------------")
	#print(tie_matrix_final)
	#print("--------------------------------------------")

	A_matrix = np.array(A_matrix, dtype = 'float').reshape(8,24)

	#for i in A_matrix:
	#	print("-----------------------------------------------------------")
	#	print(i)

	#------------------------------------------------------------------------------------------------------------------
	#-------------------------------------------Mounting A matrix: Ground control points part -------------------------
	#------------------------------------------------------------------------------------------------------------------

	zeros_matrix = np.zeros((8,6), dtype = 'float')
	GCP_matrix_coefs = np.array(GCP_list, dtype = 'float').reshape(8,6)

	GCP_matrix_final = np.concatenate((GCP_matrix_coefs, zeros_matrix, zeros_matrix, zeros_matrix), axis=1)

	#print("--------------------------------------------")
	#print(GCP_list)
	#print("--------------------------------------------")
	#print(GCP_matrix_final)

	A_matrix_final = np.concatenate((A_matrix, GCP_matrix_final), axis = 0)

	'''for i in A_matrix_final:
		print(i)
		print("--------------------------------------------")
	

	print(A_matrix_final.shape)'''


	return A_matrix_final

#------------------------------------------------------------------------------------------------------------------
#-------------------------------------------Mouting A matrix: Photo 2 -------------------------------------------
#------------------------------------------------------------------------------------------------------------------

def calcula_matA_foto2(parametros):

	a1 = diff(x,X0)
	a2 = diff(x,Y0)
	a3 = diff(x,Z0)
	a4 = diff(x,om)
	a5 = diff(x,fi)
	a6 = diff(x,kapa)
	a7 = diff(y,X0)
	a8 = diff(y,Y0)
	a9 = diff(y,Z0)
	a10 = diff(y,om)
	a11 = diff(y,fi)
	a12 = diff(y,kapa)
	a13 = diff(x, X) #tie
	a14 = diff(x, Y) #tie
	a15 = diff(x, Z) #tie
	a16 = diff(y, X) #tie
	a17 = diff(y, Y) #tie
	a18 = diff(y, Z) #tie

	coef = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18]

	cont = 0 #contador para as linhas do arquivo lido
	cont2 = 0 #Contador para controlar o calculo da lista coef
	val_a = []

	for i in parametros:
		for j in coef:
			if (cont <= 3):
				num_a = j.evalf(subs={c: 50, X0: i[10], Y0: i[11], Z0: i[12], X: i[4], Y: i[5], Z: i[6], om: i[16], fi: i[17], kapa: i[18]})
				val_a.append(num_a)
			else: #Quando for GCP, não precisa calcular derivadas em X,Y,Z
				if cont2 > 11:
					break
				else:
					num_a = j.evalf(subs={c: 50, X0: i[10], Y0: i[11], Z0: i[12], X: i[4], Y: i[5], Z: i[6], om: i[16], fi: i[17], kapa: i[18]})
					val_a.append(num_a)
					cont2 = cont2 + 1
		cont2 = 0
		cont = cont + 1
	
	tie_list = val_a[0:72]
	GCP_list = val_a[72:120]

	#------------------------------------------------------------------------------------------------------------------
	#-------------------------------------------Mouting A matrix: Tie points part ------------------------------------
	#------------------------------------------------------------------------------------------------------------------

	tie_matrix = np.array(tie_list).reshape(12,6) #
	tie_matrix_part_list = []
	A_matrix = []
	p1 = np.zeros((2,3), dtype = 'float')
	cont = 0
	cont2 = 0 #controlador da parte dos tie

	tie_matrix_final = []

	for i in tie_matrix:
		if cont != 2:
			tie_matrix_part_list.append(i)
			cont = cont + 1
		else:
			if cont2 == 0:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((p1, p1, tie_matrix_part), axis = 1) #Este cara é POE
				p2 = np.array(i).reshape(2,3) #Este cara é um diferecial em X,Y,Z
				part = np.concatenate((p2, p1, p1,p1), axis =1) #funciona adequadamente
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 == 1:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((p1, p1, tie_matrix_part), axis =1) #POE
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p2, p1, p1), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 ==2:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((p1, p1, tie_matrix_part), axis =1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p1, p2, p1), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
			elif cont2 == 3:
				tie_matrix_part = np.array(tie_matrix_part_list).reshape(2,6)
				tie_matrix_part = np.concatenate((p1, p1, tie_matrix_part), axis =1)
				p2 = np.array(i).reshape(2,3)
				part = np.concatenate((p1, p1, p1, p2), axis =1)
				tie_matrix_final = np.concatenate((tie_matrix_part, part), axis =1)
				cont2 = cont2 + 1
				cont = 0
				for k in tie_matrix_final:
					A_matrix.append(k)
				tie_matrix_part_list = []
		
		

	#print(p1)
	#print("--------------------------------------------")
	#print(tie_matrix_part)
	#print("--------------------------------------------")
	#print(p2)
	#print("--------------------------------------------")
	#print(part)
	#print("--------------------------------------------")
	#print(tie_matrix_final)
	#print("--------------------------------------------")

	A_matrix = np.array(A_matrix, dtype = 'float').reshape(8,24)

	#for i in A_matrix:
	#	print("-----------------------------------------------------------")
	#	print(i)

	#------------------------------------------------------------------------------------------------------------------
	#-------------------------------------------Mounting A matrix: Ground control points part -------------------------
	#------------------------------------------------------------------------------------------------------------------

	zeros_matrix = np.zeros((8,6), dtype = 'float')
	GCP_matrix_coefs = np.array(GCP_list, dtype = 'float').reshape(8,6)

	GCP_matrix_final = np.concatenate((zeros_matrix, GCP_matrix_coefs, zeros_matrix, zeros_matrix), axis=1)

	#print("--------------------------------------------")
	#print(GCP_list)
	#print("--------------------------------------------")
	#print(GCP_matrix_final)

	A_matrix_final_2 = np.concatenate((A_matrix, GCP_matrix_final), axis = 0)

	'''for i in A_matrix_final_2:
		print(i)
		print("--------------------------------------------")
	

	print(A_matrix_final_2.shape)'''


	return A_matrix_final_2



print('---------------------------------------------------------------------------------------------')
print('Bundle Adjustment')
print("Desenvolvido por: Natália Carvalho de Amorim")
print('Prof. Orientador: E. A. Mitishita')
print('Programa de Pós-Graduação em Ciências Geodésicas')
print('Universidade Federal do Paraná - UFPR')
print('----------------------------------------------------------------------------------------------')

#---------------------------------- ---------------------------------------------------------------
#---------------------------------- Reading the file with the initial aproximations----------------
#----------------------------------------------------- --------------------------------------------


parametros = readArqs()
parametros_copy = parametros

#---------------------------------- --------------------------------------------------------------
#---------------------------------- Converting coordinates and eliminating lens distortion--------
#----------------------------------------- and extracting L_obs------------------------------------

Lo_foto1 = []
Lo_foto2 = []

cont = 0
for i in parametros_copy:
	for j in range(len(i)):
		if (j%2 ==0) and (j <= 3):
			#print(i[j], i[j+1])
			conv = convertCoordinates(i[j], i[j+1])
			conv_corrigido = calculateRadialDist(conv[0], conv[1])
			#print(conv)
			#print(conv_corrigido)
			i[j] = conv_corrigido[0] - x0
			i[j+1] = conv_corrigido[1] - y0
			if j ==0:
				Lo_foto1.append(conv_corrigido[0] - x0)
				Lo_foto1.append(conv_corrigido[1] - y0)
			else:
				Lo_foto2.append(conv_corrigido[0]-x0)
				Lo_foto2.append(conv_corrigido[1]-y0)
			#print("----------------------------------------------------")
		elif (j > 3):
			break

Lo_foto1 = np.array(Lo_foto1, dtype = 'float').reshape(16,1)
Lo_foto2 = np.array(Lo_foto2, dtype = 'float').reshape(16,1)
L_obs = np.concatenate((Lo_foto1, Lo_foto2), axis = 0) #Observations vector



stop = False
threshold = 0.0000001
it = 1

while not stop:

	if it != 1:
		X_ant = Xv

	print('-------------------------------------------------------------------------------------------------')
	print("Iteração: ", it)
	print('-------------------------------------------------------------------------------------------------')

	#---------------------------------- Mounting the L_calc vector --------------------------------------------
	#---------------------------------- L = obs - calc -----------------------------------------------------
	#----------------------------------------------------- -------------------------------------------------

	L_calc = []
	#---------------------------------------- To photo 1 ------------------------------------------------
	for i in parametros_copy:
		L_calc.append(x.evalf(subs={c: 50, X0: i[7], Y0: i[8], Z0: i[9], X: i[4], Y: i[5], Z: i[6], om: i[13], fi: i[14], kapa: i[15]}))
		L_calc.append(y.evalf(subs={c: 50, X0: i[7], Y0: i[8], Z0: i[9], X: i[4], Y: i[5], Z: i[6], om: i[13], fi: i[14], kapa: i[15]}))

	#---------------------------------------- To photo 2 ------------------------------------------------

	for i in parametros_copy:
		L_calc.append(x.evalf(subs={c: 50, X0: i[10], Y0: i[11], Z0: i[12], X: i[4], Y: i[5], Z: i[6],  om: i[16], fi: i[17], kapa: i[18]}))
		L_calc.append(y.evalf(subs={c: 50, X0: i[10], Y0: i[11], Z0: i[12], X: i[4], Y: i[5], Z: i[6],  om: i[16], fi: i[17], kapa: i[18]}))

	L_calc = np.array(L_calc, dtype = 'float').reshape(32,1)

	L = L_obs - L_calc
	#print(L.shape)


	#---------------------------------- ---------------------------------------------------------------
	#----------------------------------Getting A matrix  --------------------------------------------
	#--------------------------------------------------------------------------------------------------


	A1 = calcula_matA_foto1(parametros_copy)
	A2 = calcula_matA_foto2(parametros_copy)

	A_final = np.concatenate((A1, A2), axis =0)
			

	#---------------------------------- ---------------------------------------------------------------
	#----------------------------------Solving the equation system ----------------------------------
	#--------------------------------------------------------------------------------------------------
	At = np.transpose(A_final)
	part1 = np.dot(At, matriz_P)
	part2 = np.dot(part1, A_final)
	N = inv(part2)

	U = np.dot(part1,L)

	'''for i in part1:
		print(i)
	print('------------------------------------------------------------------')
	for i in U:
		print(i)
	print('------------------------------------------------------------------')

	for i in N:
		print(i)'''
	Xv = np.dot(N,U)
	#print(Xv.shape)

	#---------------------------------- ---------------------------------------------------------------
	#----------------------------------Calculating residuals -------------------------------------------
	#--------------------------------------------------------------------------------------------------
	L = L_calc - L_obs
	V = np.dot(A_final, Xv) + L
	#print(V)
	#print(V.shape)

	#---------------------------------- ------------------------------------------------------------------
	#----------------------------------Verificando a diferença -------------------------------------------
	#-----------------------------------------------------------------------------------------------------

	if it !=1:
		X_dif = Xv- X_ant

	#---------------------------------- ---------------------------------------------------------------
	#----------------------------------Incrementing the initial aproximations with X values -----------
	#--------------------------------------------------------------------------------------------------
	cont = 0

	for i in parametros_copy:
		for j in range(len(i)): #vai controlar as posições de cada linha contida em parametros_copy
			if (j >=4) and (j <=6) and (cont < 4):
				if cont ==0:
					i[j] = i[j] + Xv[j+8][0]
				elif cont == 1:
					i[j] = i[j] + Xv[j+11][0]
				elif cont ==2:
					i[j] = i[j] + Xv[j+14][0]
				elif cont ==3:
					i[j] = i[j] + Xv[j+17][0]
			elif (j >=7) and (j <=9):
				i[j] = i[j] + Xv[j-7][0]
			elif (j >=10) and (j<=12):
				i[j] = i[j] + Xv[j-4][0]
			elif (j >= 13) and (j <=15):
				i[j] = i[j] + Xv[j-10][0]
			elif (j >= 16) and (j <= 18):
				i[j] = i[j] + Xv[j-7][0]
		cont = cont+1

	print(Xv)
	
	if it != 1: #Stop criteria
		if ((np.amax(Xv)) <= threshold) or (it ==10):
			stop = True
	it = it+1

	


# Print the final results

print('---------------------------Resultado Final ----------------------------------')
linha = 0
for i in parametros_copy:
	if (linha>=0) and (linha <=3):
		if linha ==0:
			print('--------------------------- Pontos Tie -----------------------------------')
			print('X: ', i[4], 'Y: ', i[5], 'Z: ', -i[6]/4)
			print('rX: ', V[12], 'rY: ', V[13], 'rZ: ', V[14], '\n')
		elif linha ==1:
			print('X: ', i[4], 'Y: ', i[5], 'Z: ', i[6]/4)
			print('rX: ', V[15], 'rY: ', V[16], 'rZ: ', V[17], '\n')
		elif linha ==2:
			print('X: ', i[4], 'Y: ', i[5], 'Z: ', i[6])
			print('rX: ', V[18], 'rY: ', V[19], 'rZ: ', V[20], '\n')
		elif linha ==3:
			print('X: ', i[4], 'Y: ', i[5], 'Z: ', i[6])
			print('rX: ', V[21], 'rY: ', V[22], 'rZ: ', V[23], '\n')
		elif linha == 4:
			print('X: ', i[4], 'Y: ', i[5], 'Z: ', i[6])
			print('rX: ', V[24], 'rY: ', V[25], 'rZ: ', V[26], '\n')
		
	elif (linha >=4):
		print('--------------------------- Foto 1 -----------------------------------')
		print('X0: ', i[7], 'Y0: ', i[8], 'Z0: ', i[9], '\nom: ', i[13], 'fi: ', i[14], 'kapa: ', i[15], '\n')
		print('rX0: ',V[0], 'rY0: ', V[1], 'rZ0: ', V[2], '\n' )
		print('--------------------------- Foto 2 -----------------------------------')
		print('X0: ', i[10], 'Y0: ', i[11], 'Z0: ', i[12],  '\nom: ', i[16], 'fi: ', i[17], 'kapa: ', i[18], '\n')
		print('rX0: ',V[6], 'rY0: ', V[7], 'rZ0: ', V[8], '\n')
		break
	linha = linha +1


#print(V)
#'\nrOm: ', V[3], 'rFi: ', V[4], 'rKapa: ', V[5], 
#'\nrOm: ', V[9], 'rFi: ', V[10], 'rKapa: ', V[11], 


