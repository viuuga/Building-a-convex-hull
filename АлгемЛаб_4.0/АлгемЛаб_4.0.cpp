/*
Алгем лабароторная работа номер 4
В этом коде реализованны все пункты
Елькин Александр М. ИВТ-11БО
Использованны дополнительные библиотеки - GLFM, GLM

ось x - красная
ось y - синяя
ось z - зелёная

Чтобы закрыть окно нажмите: esc
клавиши передвижения: w - вперёд, s - назад, a - влево, d - вправо
движение мышки - поворот камеры
колёсико на мышке - изменение градуса обзора.
клавиша E - телепортирует к первому многограннику след. нажатие ко второму многограннику, затем опять к первому и т.д.
первоначальное положение камеры - какаята из вершин 1-го многогранника
*/

#include<vector>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<locale>
#include<set>
#include <iostream>
#include <glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
struct point {
	double x, y, z;

	point operator + (point B)
	{
		point res;
		res.x = x + B.x;
		res.y = y + B.y;
		res.z = z + B.z;
		return res;
	}
	point operator - (point B)
	{
		point res;
		res.x = x - B.x;
		res.y = y - B.y;
		res.z = z - B.z;
		return res;
	}
	point operator * (double k)
	{
		point res;
		res.x = x * k;
		res.y = y * k;
		res.z = z * k;
		return res;
	}
};

const unsigned int SCR_WIDTH = 1900;
const unsigned int SCR_HEIGHT = 1000;

glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

float deltaTime = 0.0f;
float lastFrame = 0.0f;

float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float fov = 45.0f;

float yaw = -90.0f;
float pitch = 0.0f;

void processInput(GLFWwindow* window, vector<vector<pair<vector<point>, vector<pair<point, point>>>>>& edges_with_edges, int& count) // Функция для считывания нажатий на клавиатуру.
{
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	const float cameraSpeed = 10.0f * deltaTime;

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		cameraPos += cameraSpeed * cameraFront;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		cameraPos -= cameraSpeed * cameraFront;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
		if (count == 0) {
			count = 1;
			cameraPos = { edges_with_edges[0][0].second[0].first.x, edges_with_edges[0][0].second[0].first.z, edges_with_edges[0][0].second[0].first.y };
		}
		else if (count == 1 && edges_with_edges[1].size() > 0) {
			count = 0;
			cameraPos = { edges_with_edges[1][0].second[0].first.x, edges_with_edges[1][0].second[0].first.z, edges_with_edges[1][0].second[0].first.y };
		}
	}
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)  //Функция для считывания движений мышки
{
	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos;

	const float sensitivity = 0.05f;
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	yaw += xoffset;
	pitch += yoffset;

	if (pitch > 89.0f)
		pitch = 89.0f;
	if (pitch < -89.0f)
		pitch = -89.0f;

	glm::vec3 front;
	front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
	front.y = sin(glm::radians(pitch));
	front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
	cameraFront = glm::normalize(front);

	lastX = xpos;
	lastY = ypos;
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) //функция для считывания движений колёсика
{
	if (fov >= 1.0f && fov <= 120.0f)
		fov -= yoffset;
	if (fov <= 1.0f)
		fov = 1.0f;
	if (fov >= 120.0f)
		fov = 120.0f;
}

int scrian(vector<vector<pair<vector<point>, vector<pair<point, point>>>>>& edges_with_edges) // Функция создающая окно и многогранник
{
	if (!glfwInit())
	{
		std::cerr << "Failed to initialize GLFW" << std::endl;
		return -1;
	}

	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "3D Camera Movement", NULL, NULL);
	if (!window)
	{
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	cameraPos = { edges_with_edges[0][0].second[0].first.x, edges_with_edges[0][0].second[0].first.z, edges_with_edges[0][0].second[0].first.y };

	glfwMakeContextCurrent(window);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetScrollCallback(window, scroll_callback);
	int count = 0;
	while (!glfwWindowShouldClose(window))  // главный цикл в котором все рисуется и обробатывается.
	{
		processInput(window, edges_with_edges, count);

		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glm::mat4 projection = glm::perspective(glm::radians(fov), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
		glLoadMatrixf(glm::value_ptr(projection));

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glfwSetCursorPosCallback(window, mouse_callback);
		glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		glLoadMatrixf(glm::value_ptr(view));

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(-10000.0f, 0.0f, 0.0f);
		glVertex3f(10000.0f, 0.0f, 0.0f);


		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, -10000.0f, 0.0f);
		glVertex3f(0.0f, 10000.0f, 0.0f);


		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, -10000.0f);
		glVertex3f(0.0f, 0.0f, 10000.0f);
		glEnd();

		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < edges_with_edges.size(); ++i) {           // рисование многогранника.
			if (i == 1) glColor3f(0.0f, 0.0f, 0.0f);
			if (i == 0) glColor3f(1.0f, 0.0f, 1.0f);
			for (int j = 0; j < edges_with_edges[i].size(); ++j) {

				for (int k = 0; k < edges_with_edges[i][j].second.size(); ++k) {


					glVertex3f(edges_with_edges[i][j].second[k].first.x, edges_with_edges[i][j].second[k].first.z, edges_with_edges[i][j].second[k].first.y);
					glVertex3f(edges_with_edges[i][j].second[k].second.x, edges_with_edges[i][j].second[k].second.z, edges_with_edges[i][j].second[k].second.y);
				}
			}
		}


		glEnd();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}




// Функции для вывода проблем
void outprob1()
{
	cout << "Ошибка при вводе координат точек\nПрограмма завершила работу из за ошибки ввода :(";
}
void outprob2()
{
	cout << "Ошибка открытия файла\nПрограмма завершила работу из-за невозможности открыть введённый фаил :(";
}
void outprob3()
{
	cout << "Ошибка: введённое множество точек не образует многогранник\nПрограмма завершила работу из-за невозможности построения многогранника :(";
}
void outprob4()
{
	cout << "Ошибка: введён некоректный режим работы\nПрограмма завершила работу из-за невозможности определиться :(";
}
void outprob5()
{
	cout << "Ошибка при вводе колличества точек\nПрограмма завершила работу из-за ошибки ввода:(";
}
void outprob6()
{
	cout << "введено мало точек \nПрограмма завершила работу:(";
}
void outprob7()
{
	cout << "Ошибка при вводе системы уравнений\nПрограмма завершила работу:(";
}
void outprob8()
{
	cout << "Введено недостаточно систем уравнений для построения многогранника\nПрограмма завершила работу:(";
}

int error(int namber) {
	if (namber == -1) { outprob1(); return 1; };
	if (namber == -2) { outprob2(); return 1; };
	if (namber == -3) { outprob3(); return 1; };
	if (namber == -4) { outprob4(); return 1; };
	if (namber == -5) { outprob5(); return 1; };
	if (namber == -8) { outprob6(); return 1; };
	if (namber == -6) { outprob7(); return 1; };
	if (namber == -7) { outprob8(); return 1; };
	return 0;
}

// Функция считает определитель матрицы 3 на 3
double opred3_3(point A, point B, point C)
{
	return C.x * (A.y * B.z - A.z * B.y) - C.y * (A.x * B.z - A.z * B.x) + C.z * (A.x * B.y - A.y * B.x);
}

// Функция проверяет рвыенство точек
bool equality(point A, point B)
{
	if (A.x == B.x && A.y == B.y && A.z == B.z)
		return 1;
	return 0;
}

// праверяет плоскости на совподение
bool equality2(point A, point B, double d1, double d2)
{
	vector<double> a{ A.x, A.y, A.z, d1 }, b{ B.x, B.y, B.z, d2 };
	double coef = 0;
	int count = 0;

	for (int i = 0; i < 4; ++i)
	{
		if (a[i] != 0 && b[i] != 0) {
			coef = a[i] / b[i];
			break;
		}
	}

	for (int i = 0; i < 4; ++i)
	{
		if (a[i] / b[i] == coef || (a[i] == 0 && b[i] == 0)) count++;
	}

	if (count == 4) return 1;
	return 0;
}

// проверяет вектора на коллинеарность
bool equality3(point& A, point& B)
{
	vector<double> a{ A.x, A.y, A.z }, b{ B.x, B.y, B.z };
	double coef = 0;
	int count = 0;

	for (int i = 0; i < 3; ++i)
	{
		if (a[i] != 0 && b[i] != 0) {
			coef = a[i] / b[i];
			break;
		}
	}

	for (int i = 0; i < 3; ++i)
	{
		if (a[i] / b[i] == coef || (a[i] == 0 && b[i] == 0)) count++;
	}
	if (coef < 0) B = B * -1;
	if (count == 3) return 1;
	return 0;
}

// проверяет плоскости на совподение (другие входные данные
bool equality4(vector<double>& a, vector<double>& b)
{
	double coef = 0;
	int count = 0;

	for (int i = 0; i < 4; ++i)
	{
		if (a[i] != 0 && b[i] != 0) {
			coef = a[i] / b[i];
			break;
		}
	}

	for (int i = 0; i < 4; ++i)
	{
		if (a[i] / b[i] == coef || (a[i] == 0 && b[i] == 0)) count++;
	}

	if (count == 4) return 0;
	return 1;
}

//Возвращает квадрат длинны вектора
double dist2(point A)
{
	return A.x * A.x + A.y * A.y + A.z * A.z;
}

//векторное произведение, возвращает длинну
double cross(point A, point B)
{
	point C;
	C = { A.y * B.z - A.z * B.y, -(A.x * B.z - A.z * B.x), A.x * B.y - A.y * B.x };
	return sqrt(dist2(C));
}

//векторное произведение, возвращает вектор
point cross2(point A, point B)
{
	point C;
	C = { A.y * B.z - A.z * B.y, -(A.x * B.z - A.z * B.x), A.x * B.y - A.y * B.x };
	return C;
}

// находит расстояние от точки до прямой
double dist3(point A, point B, point C) {
	double h;
	h = abs(cross(A - B, A - C)) / sqrt(dist2(A - B));
	return h;
}

// находит расстояние от точки до плоскости
double dist4(point A, point B, point C, point D)
{
	point AB = B - A, AC = C - A, AD = D - A;

	double S = abs(opred3_3(AB, AC, AD));
	return S / abs(cross(AB, AC));
}

// Проверяет вектор точек на совпадения и удаляет их
void overlap(vector<point>& CH)
{
	if (CH.size() <= 1) return;
	for (int i = 0; i < CH.size() - 1; ++i)
		for (int j = i + 1; j < CH.size(); ++j) {
			if (i < CH.size() && j < CH.size()) {
				if (equality(CH[i], CH[j])) {
					CH.erase(CH.begin() + j);
					j--;
				}
			}
		}
}

//строит самый большой тетраэдр из возможных  точек
vector<point> first_CH(vector<point>& set_points)
{
	int del;
	point A, B;
	vector<point> convex_hull;

	A = set_points[0];
	// 1) берём произвольную точку 
	// 2) находим самую отдалённую от неё точку и запоминаем
	double max = 0;
	for (int i = 1; i < set_points.size(); ++i) {
		double dist = dist2(A - set_points[i]);
		if (dist > max) {
			max = dist;
			B = set_points[i];
			del = i;
		}
	}
	convex_hull.push_back(B);
	set_points.erase(set_points.begin() + del);

	// 3) первую точку удаляем так как она произвольная и находим самую отдалённую от второй, она и будет первой
	max = 0;
	for (int i = 0; i < set_points.size(); ++i) {
		double dist = dist2(B - set_points[i]);
		if (dist > max) {
			max = dist;
			A = set_points[i];
			del = i;
		}
	}
	convex_hull.push_back(A);
	set_points.erase(set_points.begin() + del);

	// 4) третья точка - самая отдалённая от прямой между первой и второй точками.
	max = 0;
	for (int i = 0; i < set_points.size(); ++i) {
		double dist = dist3(convex_hull[0], convex_hull[1], set_points[i]);
		if (dist > max) {
			max = dist;
			A = set_points[i];
			del = i;
		}
	}
	convex_hull.push_back(A);
	set_points.erase(set_points.begin() + del);

	// 5) четвёртая точка - самая отдаленная от плоскости (1ой, 2ой и 3ей)
	max = 0;
	for (int i = 0; i < set_points.size(); ++i) {
		double dist = dist4(convex_hull[0], convex_hull[1], convex_hull[2], set_points[i]);
		if (dist > max) {
			max = dist;
			A = set_points[i];
			del = i;
		}
	}
	convex_hull.push_back(A);
	set_points.erase(set_points.begin() + del);

	overlap(convex_hull);
	return convex_hull;
}

// функция в которой вводятся данные (множество точек)
int input(vector<point>& set_points)
{
	for (int i = 0; i < set_points.size(); ++i) {

		cin >> set_points[i].x >> set_points[i].y >> set_points[i].z;
		if (cin.fail()) {
			return -1;
		}
	}
	return 0;
}

// подставляет точку в уравнение плоскости и возвращает обратный знак
int plane_equation(point N, point P, int d)
{
	if (N.x * P.x + N.y * P.y + N.z * P.z + d < 0)
		return 1;
	else if (N.x * P.x + N.y * P.y + N.z * P.z + d == 0)
		return 0;
	else return -1;
}

// подставляет точку в уравнение плоскости и возвращает знак
int plane_equation_2(point N, point P, int d)
{
	if (N.x * P.x + N.y * P.y + N.z * P.z + d < 0)
		return -1;
	else if (N.x * P.x + N.y * P.y + N.z * P.z + d == 0)
		return 0;
	else return 1;
}


// отсеивает точки лежащие внутри многогранника от наружних
vector<point> out_p(vector<point>& set_points, vector<vector<point>>& edges)
{
	vector<point> out;
	// идея: если точка будет лежать одновременно под всеми гранями, то она в многограннике, и про неё можно забыть.
	// проходим по всем граням
	for (int i = 0; i < edges.size(); ++i)
	{
		double a, k, l, d;
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d = -edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z;
		for (k = 0; k < edges.size(); ++k) if (k != i) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[i].size(); ++j) {
				if (equality(edges[k][l], edges[i][j]))
					count++;
			}
			if (count == 0) break;
		}
		// строим вектор нормали так, чтобы он гарантированно указывал в наружнюю часть многогранника.
		a = plane_equation(N, edges[k][l], d);
		// подсталяем точки в уравнение плоскости
		for (int p = 0; p < set_points.size(); ++p) {
			if (plane_equation_2(N, set_points[p], d) * a > 0)
				out.push_back(set_points[p]);
		}
	}
	overlap(out);
	return out;
}

// функция выберает первую грань над которой есть хоть одна точка, потом выбирает самую дальнюю точку над гранью.
point search_pointForNewVortex(vector<point>& set_points, vector<vector<point>>& edges)
{
	point vpoint;
	int numPl;
	vector<point> abovePlan;
	for (numPl = 0; numPl < edges.size(); ++numPl)
	{
		point N = cross2(edges[numPl][0] - edges[numPl][1], edges[numPl][0] - edges[numPl][2]);
		int a, k, l, d;
		d = -edges[numPl][0].x * N.x - edges[numPl][0].y * N.y - edges[numPl][0].z * N.z;
		for (k = 0; k < edges.size(); ++k) if (k != numPl) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[numPl].size(); ++j) {
				if (equality(edges[k][l], edges[numPl][j]))
					count++;
			}
			if (count == 0) break;
		}

		a = plane_equation(N, edges[k][l], d);

		for (int p = 0; p < set_points.size(); ++p) {
			if (plane_equation_2(N, set_points[p], d) * a > 0)
				abovePlan.push_back(set_points[p]);
		}


		double max = 0;
		for (int sh = 0; sh < abovePlan.size(); ++sh) {
			double dist = dist4(edges[numPl][0], edges[numPl][1], edges[numPl][2], abovePlan[sh]);
			if (max < dist) {
				max = dist;
				vpoint = abovePlan[sh];
			}
		}
		if (abovePlan.size() > 0) break;
	}

	for (int i = 0; i < set_points.size(); ++i)
		if (equality(set_points[i], vpoint)) {
			set_points.erase(set_points.begin() + i);
			break;
		}
	return vpoint;
}

/*
разделяет грани на два случая:
1) edges - грани которые останутся без изменений.
2) edges_with_newVortex - грани над которыми новая вершина (они впоследствии будут удалены)
*/
void separation_edges(vector<vector<point>>& edges, vector<vector<point>>& edges_with_newVortex, point vortex_point)
{
	for (int i = edges.size() - 1; i >= 0; --i)
	{
		int a, k, l, d;
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d = -edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z;
		for (k = 0; k < edges.size(); ++k) if (k != i) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[i].size(); ++j) {
				if (equality(edges[k][l], edges[i][j]))
					count++;
			}
			if (count == 0) break;
		}

		a = plane_equation(N, edges[k][l], d);

		if (plane_equation_2(N, vortex_point, d) * a > 0) {
			edges_with_newVortex.push_back(edges[i]);
			edges.erase(edges.begin() + i);
		}
	}
}

//достаёт точки из граней для дальнейшей работы с ними.
void points_from_edges(vector<vector<point>>& edges, vector<vector<point>>& edges_with_newVortex, vector<point>& pt_from_egiesfp, vector<point>& point_from_edges)
{
	for (int i = 0; i < edges.size(); ++i)
		for (int j = 0; j < edges[i].size(); ++j)
			point_from_edges.push_back(edges[i][j]);

	for (int i = 0; i < edges_with_newVortex.size(); ++i)
		for (int j = 0; j < edges_with_newVortex[i].size(); ++j)
			pt_from_egiesfp.push_back(edges_with_newVortex[i][j]);
	overlap(point_from_edges);
	overlap(pt_from_egiesfp);

	for (int i = point_from_edges.size() - 1; i >= 0; --i) {
		for (int j = 0; j < pt_from_egiesfp.size(); ++j) {
			if (equality(point_from_edges[i], pt_from_egiesfp[j])) {
				point_from_edges.erase(point_from_edges.begin() + i);
				break;
			}
		}
	}
}

// строит новые грани с новой вершиной
void final_egies(vector<vector<point>>& edges, vector<point>& point_from_edges_with_newVortex, vector<point>& convex_hull, point vortex_point, vector<vector<point>>& edges_with_newVortex)
{
	int l = edges.size();
	for (int i = 0; i < l; ++i) {
		for (int i1 = 0; i1 < edges_with_newVortex.size(); ++i1) {

			int overlap_point = 0;
			vector<point> new_egie;

			for (int j = 0; j < edges[i].size(); ++j) {
				for (int j1 = 0; j1 < edges_with_newVortex[i1].size(); ++j1) {
					if (equality(edges[i][j], edges_with_newVortex[i1][j1])) {
						overlap_point++;
						new_egie.push_back(edges[i][j]);
						convex_hull.push_back(edges[i][j]);
						break;
					}
					if (overlap_point == 2) break;
				}
				if (overlap_point == 2) break;
			}

			if (overlap_point == 2) {
				new_egie.push_back(vortex_point);
				edges.push_back(new_egie);
			}
		}
	}
}

// находит грани с одинаковым уравнением плоскости и обьеденяет их.
void prov_egies(vector<vector<point>>& edges)
{

	for (int i = 0; i < edges.size() - 1; ++i) {
		int count = 0, d1, d2;
		if (i > edges.size()) break;
		point N1 = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d1 = -edges[i][0].x * N1.x - edges[i][0].y * N1.y - edges[i][0].z * N1.z;
		for (int j = i + 1; j < edges.size(); ++j) {
			if (i < edges.size() && j < edges.size()) {
				point N2 = cross2(edges[j][0] - edges[j][1], edges[j][0] - edges[j][2]);
				d2 = -edges[j][0].x * N2.x - edges[j][0].y * N2.y - edges[j][0].z * N2.z;
				if (equality2(N1, N2, d1, d2)) {
					for (int i1 = 0; i1 < edges[j].size(); ++i1) {
						edges[i].push_back(edges[j][i1]);
					}
					edges.erase(edges.begin() + j);
					overlap(edges[i]);
				}
			}
		}
	}
}

//решает систему уравнений 3 на 4 методом крамера.
bool metod_Kramera(vector<double> s1, vector<double> s2, vector<double> s3, point& res)
{
	double determinant = s1[0] * (s2[1] * s3[2] - s2[2] * s3[1]) - s1[1] * (s2[0] * s3[2] - s2[2] * s3[0]) + s1[2] * (s2[0] * s3[1] - s2[1] * s3[0]);
	if (determinant == 0) return 0;
	res.x = (s1[3] * (s2[1] * s3[2] - s2[2] * s3[1]) - s1[1] * (s2[3] * s3[2] - s2[2] * s3[3]) + s1[2] * (s2[3] * s3[1] - s2[1] * s3[3])) / determinant;
	res.y = (s1[0] * (s2[3] * s3[2] - s2[2] * s3[3]) - s1[3] * (s2[0] * s3[2] - s2[2] * s3[0]) + s1[2] * (s2[0] * s3[3] - s2[3] * s3[0])) / determinant;
	res.z = (s1[0] * (s2[1] * s3[3] - s2[3] * s3[1]) - s1[1] * (s2[0] * s3[3] - s2[3] * s3[0]) + s1[3] * (s2[0] * s3[1] - s2[1] * s3[0])) / determinant;
	return 1;
}

// вывод выпуклой оболочки/вершин многогранника.
void print_conhull(vector<point> a)
{
	char c = 65;
	int k = 0;

	cout << "Колличество вершин: " << a.size() << endl;
	for (int i = 0; i < a.size(); i++) {

		cout << c << k << ":  " << a[i].x << " " << a[i].y << " " << a[i].z << endl;
		c++;

		if (c == 91) {
			c = 65;
			k++;
		}
	}
	cout << "\n\n";
}


bool in_vec(point a, vector<point> k) {
	for (int i = 0; i < k.size(); i++) {
		if (equality(a, k[i])) {
			return 1;
		}
	}
	return 0;
}

// создаёт плоскость перпендикулярную грани и находит рёбра.
bool creat_equat_plos_and_prov(point A, point n1, point n2, vector<point>& set_points) {
	point N = cross2(n1, n2);
	double d = -A.x * N.x - A.y * N.y - A.z * N.z;
	set<int> sov;

	for (int i = 0; i < set_points.size(); ++i) {
		sov.insert(plane_equation_2(N, set_points[i], d));
	}
	if (sov.size() == 3) return 0;

	return 1;
}

// получает на вход систему уравнений и строоит выпуклую оболочку/многогранник.
int vortex_conv_hull_modeH(vector<vector<double>>& eq_sistem, vector<point>& convex_hull)
{
	int sist_size = eq_sistem.size();
	if (sist_size < 4) return -7;

	for (int s1 = 0; s1 < sist_size - 2; ++s1) {
		for (int s2 = s1 + 1; s2 < sist_size - 1; ++s2) {
			for (int s3 = s2 + 1; s3 < sist_size; ++s3) {
				point push_b;
				if (metod_Kramera(eq_sistem[s1], eq_sistem[s2], eq_sistem[s3], push_b)) {
					convex_hull.push_back(push_b);
				}
			}
		}
	}

	overlap(convex_hull);
	int count;
	point pr;
	for (int i = convex_hull.size() - 1; i >= 0; --i) {

		pr = convex_hull[i];
		count = 0;
		for (int j = 0; j < eq_sistem.size(); ++j) {
			if (pr.x * eq_sistem[j][0] + pr.y * eq_sistem[j][1] + pr.z * eq_sistem[j][2] > eq_sistem[j][3]) {
				count = 1;
				break;
			}
		}
		if (count) {
			convex_hull.erase(convex_hull.begin() + i);
		}
	}
	return 0;
}

// получает множество произвольных точек и строит выпуклую оболочку.
int quikHull(vector<point>& set_points, vector<point>& convex_hull, vector<vector<point>>& edges)
{
	// если в множестве меньше 4 точек то выходим из алгоритма
	if (set_points.size() < 4) {
		return -8;
	}
	convex_hull = first_CH(set_points);

	// если все точки лежат на одной плоскости то выходим.
	if (convex_hull.size() < 4) return -3;

	// строим начальное множество граней.
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j) {
			if (i != j) {
				edges[i].push_back(convex_hull[j]);
			}
		}

	// главный цикл
	while (set_points.size() > 0)// пока есть точки для построения многогранника.
	{
		vector<point> out; // вектор точек вне многогранника
		out = out_p(set_points, edges);
		set_points = out;
		overlap(set_points);
		if (set_points.size() == 0) break; //если не остаётся точек для расширения, то выходим.

		point vortex_point;
		vortex_point = search_pointForNewVortex(set_points, edges); // находим новую вершину.

		vector<vector<point>> edges_with_newVortex;
		separation_edges(edges, edges_with_newVortex, vortex_point);// разделяем грани на (с новой вершиной/без новой вершиныю


		vector<point> point_from_egies, point_from_edges_with_newVortex;
		points_from_edges(edges, edges_with_newVortex, point_from_edges_with_newVortex, point_from_egies);// вытаскиваем точки из граней.

		convex_hull = point_from_egies;
		final_egies(edges, point_from_edges_with_newVortex, convex_hull, vortex_point, edges_with_newVortex);// создаем грани с новой точкойю

		convex_hull.push_back(vortex_point);
		overlap(convex_hull);
		prov_egies(edges);//обьеденяем грани с одинаковым уравнением плоскости.
	}
	return 0;
}

//выводит систему уравнений задающую многогранник/выпуклую оболочку.
void equation_sistem(vector<vector<point>> edges)
{
	for (int i = 0; i < edges.size(); ++i)
	{
		double d;
		int l, a, k;
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d = (-edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z);
		for (k = 0; k < edges.size(); ++k) if (k != i) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[i].size(); ++j) {
				if (equality(edges[k][l], edges[i][j]))
					count++;
			}
			if (count == 0) break;
		}
		a = plane_equation_2(N, edges[k][l], d);
		N = N * a * -1;
		d *= a;


		if (N.x > 0) {
			cout << " + " << N.x << "x";
		}
		if (N.x < 0) {
			cout << " - " << abs(N.x) << "x";
		}

		if (N.y > 0) {
			cout << " + " << N.y << "y";
		}
		if (N.y < 0) {
			cout << " - " << abs(N.y) << "y";
		}

		if (N.z > 0) {
			cout << " + " << N.z << "z";
		}
		if (N.z < 0) {
			cout << " - " << abs(N.z) << "z";
		}
		cout << " <= " << d << endl;
	}


}

//функция выводящая граф.
vector<pair<point, vector<point>>> graph_pol(vector<vector<point>>& edges, vector<point> convex_hull, vector<pair<vector<point>, vector<pair<point, point>>>>& edges_with_edges)
{
	vector<pair<point, vector<point>>> graph_smegn(convex_hull.size());
	point vec2;

	for (int i = 0; i < convex_hull.size(); ++i) {
		graph_smegn[i].first = convex_hull[i];
	}
	for (int i = 0; i < edges.size(); ++i) {
		edges_with_edges[i].first = edges[i];
	}

	//ищем соеденённые вершины + записываем рёбра в отдельный вектор.
	for (int i = 0; i < edges.size(); ++i) {
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		double d = -edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z;

		for (int j = 0; j < edges[i].size(); ++j) {
			for (int k = 0; k < edges[i].size(); ++k) {
				if (k != j)
				{
					int znach; pair<point, point> edge;
					for (znach = 0; znach < convex_hull.size(); ++znach) if (equality(convex_hull[znach], edges[i][j])) break;


					vec2 = edges[i][j] - edges[i][k];
					if (creat_equat_plos_and_prov(edges[i][j], N, vec2, edges[i]))
					{
						edge = { edges[i][j] , edges[i][k] };
						graph_smegn[znach].second.push_back(edges[i][k]);
						edges_with_edges[i].second.push_back(edge);
					}
				}
			}
		}
	}

	//убираем совподающие рёбра по типу(AB BA)
	for (int i = 0; i < edges.size(); ++i)
	{
		for (int j = edges_with_edges[i].second.size() - 1; j >= 0; --j)
		{
			for (int k = 0; k < edges_with_edges[i].second.size(); ++k) {
				if (equality(edges_with_edges[i].second[k].first, edges_with_edges[i].second[j].second) && equality(edges_with_edges[i].second[j].first, edges_with_edges[i].second[k].second)) {
					edges_with_edges[i].second.erase(edges_with_edges[i].second.begin() + j);
					break;
				}
			}
		}
	}

	for (int i = 0; i < convex_hull.size(); ++i) {
		overlap(graph_smegn[i].second);
	}
	return graph_smegn;
}

// 1й тип вывода графа.
void output_graph1(vector<pair<point, vector<point>>>& g)
{
	cout << "\n\n";
	int k = 0;
	char c = 65;
	cout << "  ";
	for (int i = 0; i < g.size(); ++i)
	{

		cout << c << k << " ";
		c++;
		if (c == 91) {
			c = 65;
			k++;
		}
	}
	cout << endl;
	c = 65;
	for (int i = 0; i < g.size(); ++i) {

		cout << c << k << " ";
		c++;
		for (int j = 0; j < g.size(); ++j) {
			if (in_vec(g[j].first, g[i].second)) {
				cout << 1 << "  ";
			}
			else cout << 0 << "  ";
		}
		cout << endl;
		if (c == 91) {
			c = 65;
			k++;
		}
	}
}

//2й тип вывода графа.
void output_graph2(vector<pair<vector<point>, vector<pair<point, point>>>>& edges_with_edges, vector<vector<point>>& edges, vector<point>& convex_hull)
{
	cout << "\n\n";
	for (int i = 0; i < edges.size(); ++i)
	{
		double d;
		int l, a, k;
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d = (-edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z);
		for (k = 0; k < edges.size(); ++k) if (k != i) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[i].size(); ++j) {
				if (equality(edges[k][l], edges[i][j]))
					count++;
			}
			if (count == 0) break;
		}
		a = plane_equation_2(N, edges[k][l], d);
		N = N * a;
		d *= a * -1;

		cout << "Грань: ";

		if (N.x > 0) {
			cout << " + " << N.x << "x";
		}
		if (N.x < 0) {
			cout << " - " << abs(N.x) << "x";
		}
		if (N.y > 0) {
			cout << " + " << N.y << "y";
		}
		if (N.y < 0) {
			cout << " - " << abs(N.y) << "y";
		}
		if (N.z > 0) {
			cout << " + " << N.z << "z";
		}
		if (N.z < 0) {
			cout << " - " << abs(N.z) << "z";
		}
		cout << " >= " << d << endl;
		cout << "Вершины: ";
		for (int j = 0; j < edges[i].size(); ++j) {
			int k = 0, l = 0;
			char c = 65;
			while (!equality(edges[i][j], convex_hull[l])) {
				l++;
				c++;
				if (c == 91) {
					c = 65;
					k++;
				}
			}
			cout << c << k << " (" << convex_hull[l].x << " " << convex_hull[l].y << " " << convex_hull[l].z << "), ";
		}
		cout << "\nРёбра: ";
		for (int j = 0; j < edges_with_edges[i].second.size(); ++j)
		{
			int k = 0, l = 0;
			char c = 65;
			while (!equality(edges_with_edges[i].second[j].second, convex_hull[l])) {
				l++;
				c++;
				if (c == 91) {
					c = 65;
					k++;
				}
			}
			int k1 = 0, l1 = 0;
			char c1 = 65;
			while (!equality(edges_with_edges[i].second[j].first, convex_hull[l1])) {
				l1++;
				c1++;
				if (c1 == 91) {
					c1 = 65;
					k1++;
				}
			}
			cout << c1 << k1 << " " << c << k << ", ";
		}
		cout << "\n\n";
	}


}

//функция для ввода из файла.
int input2(string name, vector<point>& ch, int k, vector<pair<vector<point>, vector<pair<point, point>>>>& rebra)
{
	char mode;
	int kol, prov;
	ifstream file;
	file.open(name);
	if (!file) {
		cout << "prob1" << endl;
		return -2;
	}

	file >> mode;
	file >> kol;
	if (file.fail()) {
		return -5;
	}

	if (mode == 'V') {

		vector<vector<point>> edges(4);
		vector<point> set_points(kol), convex_hull;

		for (int i = 0; i < kol; ++i) {
			file >> set_points[i].x >> set_points[i].y >> set_points[i].z;

			if (file.fail()) {
				return -1;
			}
		}

		if (k) cout << "\n\nМножество точек 2:" << endl;
		else cout << "Множество точек 1:" << endl;
		for (int i = 0; i < kol; ++i) {
			cout << set_points[i].x << " " << set_points[i].y << " " << set_points[i].z << endl;
		}

		if (k) cout << "\n\nСистема уравнений мнгогранника 2:" << endl;
		else cout << "Система уравнений мнгогранника 1:" << endl;
		prov = quikHull(set_points, convex_hull, edges);
		if (prov != 0) return prov;
		ch = convex_hull;

		equation_sistem(edges);

		vector<pair<vector<point>, vector<pair<point, point>>>> edges_with_edges(edges.size());
		vector<pair<point, vector<point>>> graph_smegn;
		graph_smegn = graph_pol(edges, convex_hull, edges_with_edges);
		output_graph1(graph_smegn);
		output_graph2(edges_with_edges, edges, convex_hull);
		rebra = edges_with_edges;
	}
	else if (mode == 'H') {

		vector<point> convex_hull, set_points;
		vector<vector<double>> eq_sistem(kol, vector<double>(4));
		vector<vector<point>> edges(4);

		for (int i = 0; i < kol; ++i) {
			file >> eq_sistem[i][0] >> eq_sistem[i][1] >> eq_sistem[i][2] >> eq_sistem[i][3];
			if (file.fail()) {
				return -6;
			}
		}
		if (k) cout << "\n\nСистема уравнений мнгогранника 2:" << endl;
		else cout << "Система уравнений мнгогранника 1:" << endl;
		for (int i = 0; i < eq_sistem.size(); ++i) {

			cout << eq_sistem[i][0] << "x " << eq_sistem[i][1] << "y " << eq_sistem[i][2] << "z <= " << eq_sistem[i][3] << endl;
		}

		if (k) cout << "\n\nВершины мнгогранника 2:" << endl;
		else cout << "Вершины мнгогранника 1:" << endl;
		prov = vortex_conv_hull_modeH(eq_sistem, convex_hull);
		if (error(prov)) {
			return prov;
		}
		ch = convex_hull;

		set_points = convex_hull;
		convex_hull = {};
		prov = quikHull(set_points, convex_hull, edges);
		if (error(prov)) {
			return prov;
		}
		print_conhull(convex_hull);


		vector<pair<vector<point>, vector<pair<point, point>>>> edges_with_edges(edges.size());
		vector<pair<point, vector<point>>> graph_smegn;
		graph_smegn = graph_pol(edges, convex_hull, edges_with_edges);
		output_graph1(graph_smegn);
		output_graph2(edges_with_edges, edges, convex_hull);
		rebra = edges_with_edges;
	}
	else {
		return -4;
	}

	return 0;
}

// ввод из консоли.
int input3(vector<vector<double>>& eq_sistem)
{
	for (int i = 0; i < eq_sistem.size(); ++i) {

		cin >> eq_sistem[i][0] >> eq_sistem[i][1] >> eq_sistem[i][2] >> eq_sistem[i][3];
		if (cin.fail()) {
			return -6;
		}
	}


}

//вычисляет разность минковского и определяет пересекаются многогранники или нет.
bool razn_Minkovskogo(vector<point>& ch1, vector<point>& ch2)
{
	vector<vector<point>> edges(4);
	vector<point> set_points, convex_hull;
	// из каждой точки первого многогранника вычиаем все точки второго многогранника
	for (int i = 0; i < ch1.size(); ++i) {
		for (int j = 0; j < ch2.size(); ++j) {
			set_points.push_back(ch1[i] - ch2[j]);
		}
	}

	// из полученного множества точек достаём выпуклую оболочку.
	quikHull(set_points, convex_hull, edges);

	set_points = { { 0,0,0 } };
	vector<point> out;
	// проверяем на нахождение точки 0 0 0 внутри выпуклой оболочки.
	for (int i = 0; i < edges.size(); ++i)
	{
		double a, k, l, d;
		point N = cross2(edges[i][0] - edges[i][1], edges[i][0] - edges[i][2]);
		d = -edges[i][0].x * N.x - edges[i][0].y * N.y - edges[i][0].z * N.z;
		for (k = 0; k < edges.size(); ++k) if (k != i) break;

		for (l = 0; l < edges[k].size(); ++l)
		{
			int count = 0;
			for (int j = 0; j < edges[i].size(); ++j) {
				if (equality(edges[k][l], edges[i][j]))
					count++;
			}
			if (count == 0) break;
		}

		a = plane_equation(N, edges[k][l], d);

		for (int p = 0; p < set_points.size(); ++p) {
			if (plane_equation_2(N, set_points[p], d) * a > 0)
				out.push_back(set_points[p]);
		}
	}

	if (out.size() == 0) return 1;
	return 0;
}

int main(int argc, char* argv[])

{
	setlocale(LC_ALL, "Russian");

	char mode;
	int kol, prov;

	if (argc == 1) {
		cin >> mode;
		cin >> kol;
		if (cin.fail()) {
			error(-5);
			return 0;
		}
		if (mode == 'V') {
			vector<vector<point>> edges(4);
			vector<point> set_points(kol), convex_hull;
			prov = input(set_points);

			if (error(prov)) {
				return 0;
			}

			cout << "Система уравнений мнгогранника 1:" << endl;
			prov = quikHull(set_points, convex_hull, edges);
			if (error(prov)) {
				return 0;
			}


			cout << endl;
			cout << convex_hull.size() << endl;
			cout << edges.size();

			
			equation_sistem(edges);

			vector<pair<vector<point>, vector<pair<point, point>>>> edges_with_edges(edges.size());
			vector<pair<point, vector<point>>> graph_smegn;
			graph_smegn = graph_pol(edges, convex_hull, edges_with_edges);
			output_graph1(graph_smegn);
			output_graph2(edges_with_edges, edges, convex_hull);
			vector<vector<pair<vector<point>, vector<pair<point, point>>>>> rebra;

			rebra.push_back(edges_with_edges);
			scrian(rebra);
		}
		else if (mode == 'H') {
			vector<point> convex_hull, set_points;
			vector<vector<double>> eq_sistem(kol, vector<double>(4));
			vector<vector<point>> edges(4);

			prov = input3(eq_sistem);
			if (error(prov)) {
				return 0;
			}

			prov = vortex_conv_hull_modeH(eq_sistem, convex_hull);
			if (error(prov)) {
				return 0;
			}

			set_points = convex_hull;
			convex_hull = {};
			prov = quikHull(set_points, convex_hull, edges);
			if (error(prov)) {
				return 0;
			}

			vector<pair<vector<point>, vector<pair<point, point>>>> edges_with_edges(edges.size());
			vector<pair<point, vector<point>>> graph_smegn;
			graph_smegn = graph_pol(edges, convex_hull, edges_with_edges);
			output_graph1(graph_smegn);
			output_graph2(edges_with_edges, edges, convex_hull);

			vector<vector<pair<vector<point>, vector<pair<point, point>>>>> rebra;

			rebra.push_back(edges_with_edges);
			scrian(rebra);

		}
		else if (mode == 'M') {
			vector<vector<point>> edges(4);
			vector<point> set_points(kol), convex_hull, buff;
			prov = input(set_points);

			if (error(prov)) {
				return 0;
			}
			buff = set_points;
			cout << "Система уравнений мнгогранника 1:" << endl;
			prov = quikHull(set_points, convex_hull, edges);
			if (error(prov)) {
				return 0;
			}
			convex_hull = set_points;

			cout << endl;
			cout << convex_hull.size() << endl;
			cout << edges.size();


			equation_sistem(edges);

			vector<pair<vector<point>, vector<pair<point, point>>>> edges_with_edges(edges.size());
			vector<pair<point, vector<point>>> graph_smegn;
			graph_smegn = graph_pol(edges, convex_hull, edges_with_edges);
			output_graph1(graph_smegn);
			output_graph2(edges_with_edges, edges, convex_hull);
			vector<vector<pair<vector<point>, vector<pair<point, point>>>>> rebra;

			rebra.push_back(edges_with_edges);
			scrian(rebra);
		}
		else {
			error(-4);
			return 0;
		}
	}
	if (argc >= 2) {
		vector<point> convex_hull;
		vector<point> convex_hull2;
		vector<vector<pair<vector<point>, vector<pair<point, point>>>>> rebra(2);

		prov = input2(argv[1], convex_hull, 0, rebra[0]);
		if (error(prov)) return 0;

		if (argc == 3) {
			prov = input2(argv[2], convex_hull2, 1, rebra[1]);
			if (error(prov)) return 0;
			if (razn_Minkovskogo(convex_hull, convex_hull2)) {
				cout << "\nобнаружена колизия многогранников.";
			}
			else cout << "\nколизии многогранников не обнаружено.";

			scrian(rebra);
		}
		else scrian(rebra);
	}


	cout << "\n\nПрограмма завершила работу";
}