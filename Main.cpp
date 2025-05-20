#include <GL/freeglut.h>
#include <iostream>
#include <vector>
#include <list>
#include <stack>
#include <cmath>
#include <chrono>
#include <string>
#include <algorithm>





enum BorderType { LEFT, RIGHT, BOTTOM, TOP };
enum BezierMode { NONE, BERNSTEIN, DECASTELJAU };
bool drawingBezier = false;


struct Pixel{
    int x, y;
    float color[3];
};

struct Polygons {
    std::vector<std::pair<int, int>> points;
    float color[3];
    std::vector<Pixel> pixelList;
};

struct SelectionBox {
    std::pair<int, int> start;
    std::pair<int, int> end;    
    float color[3];
};


struct LinkedPoint {
    Pixel pixel;
    LinkedPoint* next;
    LinkedPoint* prev;

    LinkedPoint(int px, int py, float r, float g, float b) {
        pixel.x = px;
        pixel.y = py;
        pixel.color[0] = r;
        pixel.color[1] = g;
        pixel.color[2] = b;
        next = nullptr;
        prev = nullptr;
    }

    void addPoint(int x, int y) {
        float color[3] = {1.0f, 0.5f, 0.0f};
        LinkedPoint* newPoint = new LinkedPoint(x, y, color[0], color[1], color[2]);
        if (!next) {
            next = newPoint;
            newPoint->prev = this;
        } else {
            LinkedPoint* temp = next;
            while (temp->next) {
                temp = temp->next;
            }
            temp->next = newPoint;
            newPoint->prev = temp;
        }
    }
};

struct BezierCurve {
    LinkedPoint* head;
    std::vector<Pixel> controlPoints;
    std::vector<Pixel> curvePoints;
    BezierMode mode;
    int steps;
    double computeTime;
    float color[3];
};

struct ListChaineCourbes {
    BezierCurve list;
    ListChaineCourbes* next;
};

std::list<BezierCurve> bezierCurves;
ListChaineCourbes* bezierHead = nullptr;
std::vector<Pixel> currentControlPoints;
BezierMode currentBezierMode = NONE;
int currentBezierSteps = 1000;
const int currentBezierStepsMin = 10;
const int currentBezierStepsMax = 10000;
float currentBezierColor[3] = {1.0f, 0.5f, 0.0f};
float selectedColor[3] = {1.0f, 1.0f, 1.0f};
float colorPoint[3] = {1.0f, 0.5f, 0.0f};
int mouseX = 0, mouseY = 0;
int windowWidth = 500, windowHeight = 500;
bool starSelectionBox = false;
bool showPixel = false;
bool drawTrail = false;
bool showLine = false;
bool drawingPolygon = false;
bool drawingSelectionBox = false;
bool deleteCurve = false;
bool makeFill = false;
bool recursif = false;
bool pile = false;
bool movingPoints = false;
bool movingSoloPoint = false;
bool scanLine = false;
bool LCA = false;
std::vector<SelectionBox> selectionBoxList;
std::pair<int, int> selectionBoxStart;
std::pair<int, int> selectionBoxEnd;
std::vector<std::pair<int, int>> trail;
std::vector<std::pair<int, int>> LinePoints;
std::vector<Polygons> polygonList;
std::vector<Pixel> PixelList;
Polygons currentPolygon;
SelectionBox currentSelBox;




bool isSameColor(const float c1[3], const float c2[3], float tol = 0.01f) {
    return (std::fabs(c1[0] - c2[0]) < tol &&
            std::fabs(c1[1] - c2[1]) < tol &&
            std::fabs(c1[2] - c2[2]) < tol);
}

/**
 * @brief Lit la couleur du pixel (x, y) à l'écran et la stocke dans outColor.
 *        outColor sera normalisé dans [0..1].
 */
void getPixelColor(int x, int y, float outColor[3]) {
    GLubyte px[3];
    glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, px);
    outColor[0] = px[0] / 255.0f;
    outColor[1] = px[1] / 255.0f;
    outColor[2] = px[2] / 255.0f;
    PixelList.push_back({x, y, {outColor[0], outColor[1], outColor[2]}});
}

/**
 * @brief Dessine un pixel unique en (x, y) avec la couleur color[3].
 *        Nécessite un glFlush ou glutSwapBuffers pour être visible.
 */
void putPixel(int x, int y, const float color[3]) {
    glColor3fv(color);
    glBegin(GL_POINTS);
        glVertex2i(x, y);
    glEnd();
    glFlush(); 

}


bool isInside(const std::pair<int,int>& P, float boundaryValue, BorderType border) {
    switch (border) {
        case LEFT:
            return P.first >= boundaryValue; // x >= xmin
        case RIGHT:
            return P.first <= boundaryValue; // x <= xmax
        case BOTTOM:
            return P.second >= boundaryValue; // y >= ymin
        case TOP:
            return P.second <= boundaryValue; // y <= ymax
    }
    return false; 
}

bool pointInPolygon(int x, int y, const std::vector<std::pair<int,int>>& polygon) {
    bool inside = false;
    int n = polygon.size();
    for (int i = 0, j = n - 1; i < n; j = i++) {
        int xi = polygon[i].first, yi = polygon[i].second;
        int xj = polygon[j].first, yj = polygon[j].second;
        if (((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (float)(yj - yi) + xi))
            inside = !inside;
    }
    return inside;
}


void floodFillRecursive(const std::pair<int,int>& P,float boundaryColor[3])
{
    int x = P.first;
    int y = P.second;

    float currentColor[3];
    getPixelColor(x, y, currentColor);

    if (isSameColor(currentColor, boundaryColor)) {
        return; 
    }

    putPixel(x, y, boundaryColor);

    floodFillRecursive(std::make_pair(x+1, y), boundaryColor);
    floodFillRecursive(std::make_pair(x-1, y), boundaryColor);
    floodFillRecursive(std::make_pair(x, y+1), boundaryColor);
    floodFillRecursive(std::make_pair(x, y-1), boundaryColor);

    glutSwapBuffers();
}

void floodFillStack(std::pair<int,int>& P, float boundaryColor[3]){
    std::stack<std::pair<int, int>> pile;
    pile.push({P.first, P.second});

    while(!pile.empty()){
        auto[x,y] = pile.top();
        pile.pop();

        float currentColor[3];
        getPixelColor(x, y, currentColor);

        if (isSameColor(currentColor, boundaryColor)) {
            continue; 
        }

        putPixel(x, y, boundaryColor);

        pile.push(std::make_pair(x+1, y));
        pile.push(std::make_pair(x-1, y));
        pile.push(std::make_pair(x, y+1));
        pile.push(std::make_pair(x, y-1));
    }

}

void scanlineFloodFill(int startX, int startY, float boundaryColor[3]) {
    std::stack<std::pair<int,int>> pile;
    pile.push({startX, startY});

    while (!pile.empty()) {
        auto top = pile.top();
        int x = top.first;
        int y = top.second;
        pile.pop();

        int xLeft = x;
        float currentColor[3];
        while (xLeft >= 0) {
            getPixelColor(xLeft, y, currentColor);
            if (isSameColor(currentColor, boundaryColor))
                break;
            xLeft--;
        }
        xLeft++; 

        bool spanUp = false;
        bool spanDown = false;

        while (xLeft < windowWidth) {
            getPixelColor(xLeft, y, currentColor);
            if (isSameColor(currentColor, boundaryColor))
                break;

            putPixel(xLeft, y, boundaryColor);

            if (y > 0) {
                float aboveColor[3];
                getPixelColor(xLeft, y - 1, aboveColor);
                if (!(isSameColor(aboveColor, boundaryColor))) {
                    if (!spanUp) {
                        pile.push({xLeft, y - 1});
                        spanUp = true;
                    }
                } else {
                    spanUp = false;
                }
            }

            if (y < windowHeight - 1) {
                float belowColor[3];
                getPixelColor(xLeft, y + 1, belowColor);
                if (!(isSameColor(belowColor, boundaryColor))) {
                    if (!spanDown) {
                        pile.push({xLeft, y + 1});
                        spanDown = true;
                    }
                } else {
                    spanDown = false;
                }
            }
            xLeft++;
        }
    }
}

std::vector<int> pascalLine(int n) {
    std::vector<int> line(n + 1, 1);
    for (int i = 1; i < n; ++i) {
        line[i] = line[i - 1] * (n - i + 1) / i;
    }
    return line;
}

// Fonction pour calculer un point de la courbe de Bézier pour un t donné
Pixel bezierPoint(const std::vector<Pixel>& controlPoints, double t) {
    int n = controlPoints.size() - 1;
    std::vector<int> C = pascalLine(n);
    Pixel result = {0, 0, {0, 0, 0}};

    for (int i = 0; i <= n; ++i) {
        double coeff = C[i] * pow(1 - t, n - i) * pow(t, i);
        result.x += coeff * controlPoints[i].x;
        result.y += coeff * controlPoints[i].y;
    }

    return result;
}

// Fonction pour générer la courbe de Bézier en utilisant des pas de t
std::vector<Pixel> bezierCurve(const LinkedPoint* head, int steps = 100) {
    std::vector<Pixel> curvePoints;
    const LinkedPoint* current = head;
    std::vector<Pixel> controlPoints;
    while (current != nullptr) {
        controlPoints.push_back(current->pixel);
        current = current->next;
    }

    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        curvePoints.push_back(bezierPoint(controlPoints, t));
    }
    return curvePoints;
}

Pixel deCasteljauIterative(const std::vector<Pixel>& controlPoints, double t) {
    std::vector<double> px, py;
    for (const Pixel& p : controlPoints) {
        px.push_back(p.x);
        py.push_back(p.y);
    }
    int n = controlPoints.size();
    for (int k = 1; k < n; ++k) {
        for (int i = 0; i < n - k; ++i) {
            px[i] = (1 - t) * px[i] + t * px[i + 1];
            py[i] = (1 - t) * py[i] + t * py[i + 1];
        }
    }
    Pixel result;
    result.x = static_cast<int>(px[0]);
    result.y = static_cast<int>(py[0]);
    result.color[0] = controlPoints[0].color[0];
    result.color[1] = controlPoints[0].color[1];
    result.color[2] = controlPoints[0].color[2];
    return result;
}

std::vector<Pixel> bezierCurveDeCasteljau(const std::vector<Pixel>& controlPoints, int steps = 100) {
    std::vector<Pixel> curvePoints;
    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        curvePoints.push_back(deCasteljauIterative(controlPoints, t));
    }
    return curvePoints;
}


std::pair<int,int> computeIntersectionSutherland(
    const std::pair<int,int>& S1, 
    const std::pair<int,int>& S2, 
    float boundaryValue, 
    BorderType border
) {
    float x1 = (float)S1.first,  y1 = (float)S1.second;
    float x2 = (float)S2.first,  y2 = (float)S2.second;
    float x, y;
    
    if (border == LEFT || border == RIGHT) {
        if (std::abs(x2 - x1) < 1e-5f) {
            return S1; 
        }
        float t = (boundaryValue - x1) / (x2 - x1);
        x = boundaryValue;
        y = y1 + t * (y2 - y1);
    } 
    else {
        if (std::abs(y2 - y1) < 1e-5f) {
            return S1;
        }
        float t = (boundaryValue - y1) / (y2 - y1);
        y = boundaryValue;
        x = x1 + t * (x2 - x1);
    }

    return std::make_pair((int)x, (int)y);
}

std::vector<std::pair<int,int>> clipWithBorderSutherland(
    const std::vector<std::pair<int,int>>& inPolygon,
    float boundaryValue,
    BorderType border
) {
    std::vector<std::pair<int,int>> outPolygon;
    for (size_t i = 0; i < inPolygon.size(); i++) {
        std::pair<int,int> current = inPolygon[i];
        std::pair<int,int> next = inPolygon[(i + 1) % inPolygon.size()];

        bool currInside = isInside(current, boundaryValue, border);
        bool nextInside = isInside(next, boundaryValue, border);

        if (currInside && nextInside) {
            outPolygon.push_back(next);
        } 
        else if (currInside && !nextInside) {
            outPolygon.push_back(computeIntersectionSutherland(current, next, boundaryValue, border));
        }
        else if (!currInside && nextInside) {
            outPolygon.push_back(computeIntersectionSutherland(current, next, boundaryValue, border));
            outPolygon.push_back(next);
        }
    }
    return outPolygon;
}

Polygons clipPolygonSutherland(const Polygons& inputPoly, const SelectionBox& box) {
    Polygons output;

    float xmin = std::min(box.start.first, box.end.first);
    float xmax = std::max(box.start.first, box.end.first);
    float ymin = std::min(box.start.second, box.end.second);
    float ymax = std::max(box.start.second, box.end.second);

    std::vector<std::pair<int,int>> clippedPoints = inputPoly.points;

    clippedPoints = clipWithBorderSutherland(clippedPoints, xmin, LEFT);

    clippedPoints = clipWithBorderSutherland(clippedPoints, xmax, RIGHT);

    clippedPoints = clipWithBorderSutherland(clippedPoints, ymin, BOTTOM);

    clippedPoints = clipWithBorderSutherland(clippedPoints, ymax, TOP);

    output.points = clippedPoints;

    output.color[0] = 1.0f;
    output.color[1] = 0.0f;
    output.color[2] = 1.0f;

    return output;
}

std::pair<int,int> computeIntersectionCyrusBeck(
    const std::pair<int,int>& S1, 
    const std::pair<int,int>& S2, 
    float boundaryValue, 
    BorderType border
) {
    float x1 = static_cast<float>(S1.first), y1 = static_cast<float>(S1.second);
    float x2 = static_cast<float>(S2.first), y2 = static_cast<float>(S2.second);
    float dx = x2 - x1, dy = y2 - y1;
    
    // Définir la normale N et un point F sur le bord en fonction du type
    float Nx = 0.0f, Ny = 0.0f, Fx = 0.0f, Fy = 0.0f;
    switch(border) {
        case LEFT:
            // Pour x >= boundaryValue, la normale extérieure est (-1, 0)
            Nx = -1.0f; Ny = 0.0f;
            Fx = boundaryValue; Fy = 0.0f; // F peut être (xmin, 0)
            break;
        case RIGHT:
            // Pour x <= boundaryValue, la normale extérieure est (1, 0)
            Nx = 1.0f; Ny = 0.0f;
            Fx = boundaryValue; Fy = 0.0f;
            break;
        case BOTTOM:
            // Pour y >= boundaryValue, la normale extérieure est (0, -1)
            Nx = 0.0f; Ny = -1.0f;
            Fx = 0.0f; Fy = boundaryValue;
            break;
        case TOP:
            // Pour y <= boundaryValue, la normale extérieure est (0, 1)
            Nx = 0.0f; Ny = 1.0f;
            Fx = 0.0f; Fy = boundaryValue;
            break;
    }
    
    // Calcul du produit scalaire D · N
    float DdotN = dx * Nx + dy * Ny;
    if (std::abs(DdotN) < 1e-5f) {
        // Le segment est pratiquement parallèle au bord : on retourne S1 par défaut
        return S1;
    }
    // Calcul de t avec t = (N · (F - S1)) / (N · D)
    float numerator = (Fx - x1) * Nx + (Fy - y1) * Ny;
    float t = numerator / DdotN;
    
    float xi = x1 + t * dx;
    float yi = y1 + t * dy;
    return std::make_pair(static_cast<int>(xi), static_cast<int>(yi));
}

// Fonction de clipping par un bord avec Cyrus-Beck
std::vector<std::pair<int,int>> clipWithBorderCyrusBeck(
    const std::vector<std::pair<int,int>>& inPolygon,
    float boundaryValue,
    BorderType border
) {
    std::vector<std::pair<int,int>> outPolygon;
    for (size_t i = 0; i < inPolygon.size(); i++) {
        std::pair<int,int> current = inPolygon[i];
        std::pair<int,int> next = inPolygon[(i + 1) % inPolygon.size()];

        bool currInside = isInside(current, boundaryValue, border);
        bool nextInside = isInside(next, boundaryValue, border);

        if (currInside && nextInside) {
            outPolygon.push_back(next);
        } 
        else if (currInside && !nextInside) {
            outPolygon.push_back(computeIntersectionCyrusBeck(current, next, boundaryValue, border));
        }
        else if (!currInside && nextInside) {
            outPolygon.push_back(computeIntersectionCyrusBeck(current, next, boundaryValue, border));
            outPolygon.push_back(next);
        }
        // Si les deux points sont à l'extérieur, aucun point n'est ajouté
    }
    return outPolygon;
}

// Fonction principale de clipping du polygone
Polygons clipPolygonCyrusBeck(const Polygons& inputPoly, const SelectionBox& box) {
    Polygons output;

    float xmin = std::min(box.start.first, box.end.first);
    float xmax = std::max(box.start.first, box.end.first);
    float ymin = std::min(box.start.second, box.end.second);
    float ymax = std::max(box.start.second, box.end.second);

    std::vector<std::pair<int,int>> clippedPoints = inputPoly.points;

    clippedPoints = clipWithBorderCyrusBeck(clippedPoints, xmin, LEFT);
    clippedPoints = clipWithBorderCyrusBeck(clippedPoints, xmax, RIGHT);
    clippedPoints = clipWithBorderCyrusBeck(clippedPoints, ymin, BOTTOM);
    clippedPoints = clipWithBorderCyrusBeck(clippedPoints, ymax, TOP);

    output.points = clippedPoints;

    // Couleur d'exemple
    output.color[0] = 1.0f;
    output.color[1] = 0.0f;
    output.color[2] = 1.0f;

    return output;
}



void colorMenuCallBack(int value){
    switch(value){
        case 1:
            selectedColor[0] = 0.0f; selectedColor[1] = 1.0f; selectedColor[2] = 0.0f; // Vert
            std::cout << "Couleur verte selectionnee" << std::endl;
            break;
        case 2:
            selectedColor[0] = 1.0f; selectedColor[1] = 0.0f; selectedColor[2] = 0.0f; // Rouge
            std::cout << "Couleur rouge selectionnee" << std::endl;
            break;
        case 3:
            selectedColor[0] = 0.0f; selectedColor[1] = 0.0f; selectedColor[2] = 1.0f; // Bleu
            std::cout << "Couleur bleue selectionnee" << std::endl;
            break;
    }
}

void FillMenuCallBack(int value){
    switch(value){
        case 1:
            makeFill = true;
            recursif = true;
            break;
        case 2:
            makeFill= true;
            pile = true;
            break;
        case 3:
            makeFill = true;
            scanLine = true;
            break;
    }
}

void FenMenuCallBack(int value){
    std::vector<Polygons> newPolygon;
    bool clippedFlag = false;
    switch(value){
        case 1:
            if (!polygonList.empty() && !selectionBoxList.empty()) {
                for(Polygons p : polygonList){
                    for(SelectionBox sb : selectionBoxList){

                        Polygons clipped = clipPolygonSutherland(p, sb);
                        if (!clipped.points.empty()) {
                            newPolygon.push_back(clipped);
                            clippedFlag = true;
                        }
                    }
                    if (!clippedFlag)
                    {
                        newPolygon.push_back(p);
                    }else{
                        clippedFlag = false;
                    }
                    
                }
            }
            polygonList.clear();
            selectionBoxList.clear();
            polygonList = newPolygon;
            break;
        case 2:
            if (!polygonList.empty() && !selectionBoxList.empty()) {
                for(Polygons p : polygonList){
                    for(SelectionBox sb : selectionBoxList){

                        Polygons clipped = clipPolygonCyrusBeck(p, sb);
                        if (!clipped.points.empty()) {
                            newPolygon.push_back(clipped);
                            clippedFlag = true;
                        }
                    }
                    if (!clippedFlag)
                    {
                        newPolygon.push_back(p);
                    }else{
                        clippedFlag = false;
                    }
                    
                }
            }
            polygonList.clear();
            selectionBoxList.clear();
            polygonList = newPolygon;
            break;
        case 3:
            selectedColor[0] = 0.0f; selectedColor[1] = 0.0f; selectedColor[2] = 1.0f; // Bleu
            std::cout << "Couleur bleue selectionnee" << std::endl;
            break;
    }
}



void menuCallBack(int value){
    switch (value) {
        case 1:
            break;
        case 2:
            showPixel = true;
            showLine = true;
            break;
        case 3:
            drawingSelectionBox = true;
            break;
        case 4:
            movingPoints = true;
            drawingBezier = false;
            break;
        case 5:
            movingPoints = false;
            drawingBezier = false;
            deleteCurve = true;
            break;

    }
}


void drawSquare(int x, int y, int size) {
    glBegin(GL_QUADS);
    glVertex2i(x + size / 2, y - size / 2);
    glVertex2i(x - size / 2, y - size / 2);
    glVertex2i(x - size / 2, y + size / 2);
    glVertex2i(x + size / 2, y + size / 2);
    glEnd();
}

/**
 * @brief Dessine une série de lignes connectées basées sur les points stockés dans LinePoints.
 * 
 * Cette fonction utilise OpenGL pour dessiner des lignes entre des points consécutifs dans le vecteur LinePoints.
 * S'il y a moins de deux points dans le vecteur, la fonction retourne sans rien dessiner.
 * 
 * @note Suppose que LinePoints est un vecteur de paires, où chaque paire représente les coordonnées (x, y) d'un point.
 */
void drawLine() {
    if (LinePoints.size() < 2) return; 
    glBegin(GL_LINES);
    for (size_t i = 0; i < LinePoints.size() - 1; ++i) {
        glVertex2i(LinePoints[i].first, LinePoints[i].second);
        glVertex2i(LinePoints[i + 1].first, LinePoints[i + 1].second);
    }
    glEnd();
}

void drawPolygon() {
    if (currentPolygon.points.size() < 3) return; 
    glBegin(GL_LINE_LOOP);
    for (const auto& point : currentPolygon.points) {
        glVertex2i(point.first, point.second);
    }
    glEnd();
}

void drawSelectionBox() {
    glBegin(GL_LINE_LOOP);
    glVertex2i(currentSelBox.start.first, currentSelBox.start.second);
    glVertex2i(currentSelBox.end.first, currentSelBox.start.second);
    glVertex2i(currentSelBox.end.first, currentSelBox.end.second);
    glVertex2i(currentSelBox.start.first, currentSelBox.end.second);
    glEnd();

}

void polygonStay(Polygons p){
    glBegin(GL_LINE_LOOP);
    glColor3fv(p.color);
    for (const auto& point : p.points) {
        glVertex2i(point.first, point.second);
    }
    if(!p.pixelList.empty()){
        for(Pixel p : p.pixelList){
            glColor3fv(p.color);
            glVertex2i(p.x, p.y);
        }
    }
    glEnd();
}

void SelBoxStay(SelectionBox sb){
    glBegin(GL_LINE_LOOP);
    glColor3fv(sb.color);
    glVertex2i(sb.start.first, sb.start.second);
    glVertex2i(sb.end.first, sb.start.second);
    glVertex2i(sb.end.first, sb.end.second);
    glVertex2i(sb.start.first, sb.end.second);
    glEnd();
}

void displayCallBack(){
    glClear(GL_COLOR_BUFFER_BIT);
    if(showPixel){
        glColor3fv(selectedColor);
        drawSquare(mouseX, mouseY, 5);
    }

    if (!trail.empty()) {
        glColor3fv(selectedColor);
        glBegin(GL_POINTS);
        for (const auto& point : trail) {
            glVertex2i(point.first, point.second);
        }
        glEnd();
    }

    if (showLine && !LinePoints.empty()) {
        glColor3fv(selectedColor);
        drawLine();
    }

    if (drawingPolygon && !currentPolygon.points.empty()) {
        glColor3fv(selectedColor);
        drawPolygon();
    }

    if (!currentControlPoints.empty()) {
        glColor3fv(colorPoint); // Par exemple orange
        glPointSize(7);
        glBegin(GL_POINTS);
        for (const Pixel& pt : currentControlPoints) {
            glVertex2i(pt.x, pt.y);
        }
        glEnd();
        glPointSize(1);
    }

    if (currentControlPoints.size() >= 2) {
        glColor3fv(selectedColor); // ou une couleur différente
        glBegin(GL_LINE_STRIP);
        for (const Pixel& pt : currentControlPoints) {
            glVertex2i(pt.x, pt.y);
        }
        glEnd();
    }


    // Affichage des courbes de Bézier
    ListChaineCourbes* current = bezierHead;
    while (current != nullptr) {
        glColor3fv(current->list.color);
        glBegin(GL_LINE_STRIP);
        for (const Pixel& pt : current->list.curvePoints) {
            glVertex2i(pt.x, pt.y);
        }
        glEnd();
        glPointSize(7);
        glBegin(GL_POINTS);
        LinkedPoint* currentPoint = current->list.head;
        while (currentPoint != nullptr) {
            glColor3fv(currentPoint->pixel.color);
            glVertex2i(currentPoint->pixel.x, currentPoint->pixel.y);
            currentPoint = currentPoint->next;
        }
        glEnd();
        current = current->next;
    }
    // for (const auto& curve : bezierCurves) {
    //     glColor3fv(curve.color);
    //     glBegin(GL_LINE_STRIP);
    //     LinkedPoint* current = curve.head;
    //     for (const Pixel& pt : curve.curvePoints) {
    //         glVertex2i(pt.x, pt.y);
    //     }
    //     glEnd();
    
    //     // Affiche aussi les points de contrôle (optionnel)
    //     glColor3fv(curve.color);
    //     glPointSize(7);
    //     glBegin(GL_POINTS);
    //     current = curve.head;
    //     while (current != nullptr) {
    //         glVertex2i(current->pixel.x, current->pixel.y);
    //         current = current->next;
    //     }
    //     glEnd();
    //     glPointSize(1);
    // }

    if(!polygonList.empty()){
        for(Polygons p : polygonList){
            polygonStay(p);
        }
    }

    if (drawingSelectionBox && starSelectionBox) {
        glColor3fv(currentSelBox.color);
        drawSelectionBox();
    }
    if(!selectionBoxList.empty()){
        for(SelectionBox sb : selectionBoxList){
            SelBoxStay(sb);
        }
    }

    glColor3f(0.0, 0.0, 0.0);
    glRasterPos2i(10, windowHeight - 40);
    std::string stepText = "Pas (steps): " + std::to_string(currentBezierSteps);
    for (char c : stepText) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c);

    glutSwapBuffers();
}

void passiveMotionCallback(int x, int y){
    if(showPixel){
        mouseX = x;
        mouseY = windowHeight - y;
    }
    if (drawingSelectionBox && starSelectionBox) {
        currentSelBox.end = std::make_pair(x, windowHeight - y);
    }
    glutPostRedisplay();
}


void motionCallBack(int x, int y){
    if(showPixel){
        mouseX = x;
        mouseY = windowHeight - y;
    }

    if (drawingSelectionBox && starSelectionBox) {
        currentSelBox.end = std::make_pair(x, windowHeight - y);
    }
    glutPostRedisplay();
}

void mouseCallBack(int input, int state, int x, int y){
    if (input == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
         if (drawingBezier) {
            Pixel clickedPixel = {x, windowHeight - y, {currentBezierColor[0], currentBezierColor[1], currentBezierColor[2]}};
            currentControlPoints.push_back(clickedPixel);
            glutPostRedisplay();
        }
        if (showLine && showPixel) {
                drawTrail = true;
                trail.emplace_back(x, windowHeight - y);
                LinePoints.emplace_back(x, windowHeight - y); 
                currentPolygon.points.emplace_back(x, windowHeight - y);
            }
        if (drawingSelectionBox) {
            if(!starSelectionBox){
                currentSelBox.color[0] = selectedColor[0];
                currentSelBox.color[1] = selectedColor[1];
                currentSelBox.color[2] = selectedColor[2];
                std::cout << "start" << std::endl;
                currentSelBox.start = std::make_pair(x, windowHeight - y);
                starSelectionBox = true;
            }else{
                currentSelBox.end = std::make_pair(x, windowHeight - y);
                selectionBoxList.push_back(currentSelBox);
                currentSelBox.start = std::make_pair(0, 0);
                currentSelBox.end = std::make_pair(0, 0);
                starSelectionBox = false;
                drawingSelectionBox = false;
            }
        }
        if (deleteCurve) {
            std::cout << "Suppression de la courbe" << std::endl;

            ListChaineCourbes* currentCurve = bezierHead;
            std::cout << "x: " << x << " y: " << (windowHeight - y) << std::endl;
            while (currentCurve != nullptr) {
                if (std::any_of(currentCurve->list.curvePoints.begin(), currentCurve->list.curvePoints.end(), 
                    [x, y](const Pixel& p) { return std::abs(p.x - x) < 5 && std::abs(p.y - (windowHeight - y)) < 5; })) {
                    // Supprimer la courbe
                    std::cout << "Courbe supprimée" << std::endl;
                    if (currentCurve == bezierHead) {
                        bezierHead = currentCurve->next;
                    } else {
                        ListChaineCourbes* prevCurve = bezierHead;
                        while (prevCurve->next != currentCurve) {
                            prevCurve = prevCurve->next;
                        }
                        prevCurve->next = currentCurve->next;
                    }
                    deleteCurve = false;
                }
                currentCurve = currentCurve->next;
            }
        }
         // Déplacement des points de contrôle
         // std::cout << "x: " << x << " y: " << y << std::endl;
         // std::cout << "mouseX: " << mouseX << " mouseY: " << mouseY << std::endl;
        if (movingSoloPoint) {
            ListChaineCourbes* currentCurve = bezierHead;
            while (currentCurve != nullptr) {
                LinkedPoint* currentPoint = currentCurve->list.head;
                while (currentPoint != nullptr) {
                    if (currentPoint->pixel.color[0] == 0.0f && currentPoint->pixel.color[1] == 1.0f && currentPoint->pixel.color[2] == 0.0f) {
                        currentPoint->pixel.x = x;
                        currentPoint->pixel.y = windowHeight - y;
                        currentPoint->pixel.color[0] = 1.0f;
                        currentPoint->pixel.color[1] = 0.5f;
                        currentPoint->pixel.color[2] = 0.0f;
                        movingSoloPoint = false;
                        auto start = std::chrono::high_resolution_clock::now();
                        if (currentCurve->list.mode == BERNSTEIN)
                        currentCurve->list.curvePoints = bezierCurve(currentCurve->list.head, currentCurve->list.steps);
                        else
                        currentCurve->list.curvePoints = bezierCurveDeCasteljau(currentCurve->list.controlPoints, currentCurve->list.steps);
                        auto end = std::chrono::high_resolution_clock::now();
                        currentCurve->list.computeTime = std::chrono::duration<double, std::milli>(end - start).count();
                        glutPostRedisplay();
                        movingPoints = true;
                        break;
                    }
                    currentPoint = currentPoint->next;
                }
                currentCurve = currentCurve->next;
            }
        }
        else if (movingPoints) {
            std::cout << "Déplacement des points de contrôle" << std::endl;
            ListChaineCourbes* currentCurve = bezierHead;
            while (currentCurve != nullptr) {
                LinkedPoint* currentPoint = currentCurve->list.head;
                while (currentPoint != nullptr) {
                    // std::cout << "x: " << x << " y: " << y << std::endl;
                    // std::cout << "currentPoint.x: " << currentPoint->pixel.x << " currentPoint.y: " << currentPoint->pixel.y << std::endl;
                    // std::cout << "Différence x: " << std::abs(currentPoint->pixel.x - x) << " Différence y: " << std::abs(currentPoint->pixel.y - (windowHeight - y)) << std::endl;
                    if (std::abs(currentPoint->pixel.x - x) <= 5 && std::abs(currentPoint->pixel.y - (windowHeight - y)) <= 5) {
                        std::cout << "Point de contrôle déplacé" << std::endl;
                        currentPoint->pixel.color[0] = 0.0f;
                        currentPoint->pixel.color[1] = 1.0f;
                        currentPoint->pixel.color[2] = 0.0f;
                        movingSoloPoint = true;
                        movingPoints = false;
                    }
                    currentPoint = currentPoint->next;
                }
                currentCurve = currentCurve->next;
                glutPostRedisplay();
            }
        }
        if (makeFill && !drawingPolygon) {
        int adjustedY = windowHeight - y; 
            if(recursif){
                for (auto& poly : polygonList) {
                    if (pointInPolygon(x, adjustedY, poly.points)) {
                        std::cout << "Remplissage du polygone à partir du point (" 
                                << x << ", " << adjustedY << ")" << std::endl;
                        floodFillRecursive(std::make_pair(x, adjustedY), poly.color);
                        poly.pixelList = PixelList;
                        PixelList.clear();
                        makeFill = false;
                        recursif = false;
                        break;
                    }
                }
            }else if(pile){
                for (auto& poly : polygonList) {
                    if (pointInPolygon(x, adjustedY, poly.points)) {
                        std::cout << "Remplissage du polygone à partir du point (" 
                                << x << ", " << adjustedY << ")" << std::endl;
                        std::pair<int, int> fillPoint = std::make_pair(x, adjustedY);
                        floodFillStack(fillPoint, poly.color);
                        poly.pixelList = PixelList;
                        PixelList.clear();
                        makeFill = false;
                        pile = false;
                        break;
                    }
                }  
            }else if(scanLine){
                for (auto& poly : polygonList) {
                    if (pointInPolygon(x, adjustedY, poly.points)) {
                        std::cout << "Remplissage du polygone à partir du point (" 
                                << x << ", " << adjustedY << ")" << std::endl;
                       scanlineFloodFill(x, adjustedY, poly.color);
                        poly.pixelList = PixelList;
                        PixelList.clear();
                        makeFill = false;
                        scanLine = false;
                        break;
                    }
                }  
            } 
        }
    } else if (input == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        drawTrail = false;
    } 
    glutPostRedisplay();
}

void reshapeCallback(int width, int height) {
    windowWidth = width;
    windowHeight = height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, 0, height, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


void keyboard(unsigned char key, int x, int y){
    switch (key) {
        case 13:
            currentPolygon.color[0] = selectedColor[0];
            currentPolygon.color[1] = selectedColor[1];
            currentPolygon.color[2] = selectedColor[2];
            polygonList.push_back(currentPolygon);
            showPixel = false;
            drawingPolygon = false;
            showLine = false;
            currentPolygon.points.clear();
            LinePoints.clear();
            trail.clear();

            break;
        case 27:
            polygonList.clear();
            selectionBoxList.clear();
            currentControlPoints.clear();
            currentBezierMode = NONE;
            drawingBezier = false;
            break;
        case '+':
        case '=': // Pour le + du clavier principal
            if (currentBezierSteps < currentBezierStepsMax)
                currentBezierSteps += 10;
            // Recalcule la courbe si déjà affichée
            std::cout << "Pas augmenté: " << currentBezierSteps << std::endl;
            glutPostRedisplay();
            break;
        case '-':
        case '_': // Pour le - du clavier principal
            if (currentBezierSteps > currentBezierStepsMin)
                currentBezierSteps -= 10;
            std::cout << "Pas diminué: " << currentBezierSteps << std::endl;
            glutPostRedisplay();
            break;
    }
    glutPostRedisplay();
}

void timerCallback(int value) {
    glutPostRedisplay(); 
    glutTimerFunc(16, timerCallback, 0);
}

int main(int argc, char **argv){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500,500);
    glutCreateWindow("Main Window");


    int colorMenu = glutCreateMenu(colorMenuCallBack);

    glutAddMenuEntry("Vert", 1);
    glutAddMenuEntry("Rouge", 2);
    glutAddMenuEntry("Bleu", 3);

    int FillMenu = glutCreateMenu(FillMenuCallBack);

    glutAddMenuEntry("Recursif", 1);
    glutAddMenuEntry("Pile", 2);
    glutAddMenuEntry("ScanLine", 3);
    glutAddMenuEntry("LCA", 4);

    int FenMenu = glutCreateMenu(FenMenuCallBack);
    glutAddMenuEntry("Sutherland - Hodgman", 1);
    glutAddMenuEntry("CyrusBeck", 2);

    int PosePointsMenu = glutCreateMenu([](int value){
        if (value == 1) {
            currentBezierMode = NONE;
            drawingBezier = true;
            currentControlPoints.clear();
            std::cout << "Mode : Pose des points de contrôle activé." << std::endl;
        }
    });
    glutAddMenuEntry("Poser points de Bézier", 1);

    // Menu pour sélectionner l'algorithme de tracé
    int TraceBezierMenu = glutCreateMenu([](int value){
        if (currentControlPoints.size() < 2) return;

        BezierCurve newCurve;
        newCurve.head = NULL;
        LinkedPoint* current = NULL;
        for (const Pixel& p : currentControlPoints) {
            if (newCurve.head == NULL) {
                newCurve.head = new LinkedPoint(p.x, p.y, p.color[0], p.color[1], p.color[2]);
                current = newCurve.head;
            } else {
                current->addPoint(p.x, p.y);
                current = current->next;
            }
        }
        newCurve.mode = (value == 1 ? BERNSTEIN : DECASTELJAU);
        newCurve.steps = currentBezierSteps;
        newCurve.color[0] = newCurve.head->pixel.color[0];
        newCurve.color[1] = newCurve.head->pixel.color[1];
        newCurve.color[2] = newCurve.head->pixel.color[2];
        // newCurve.color[0] = currentBezierColor[0];
        // newCurve.color[1] = currentBezierColor[1];
        // newCurve.color[2] = currentBezierColor[2];
    
        auto start = std::chrono::high_resolution_clock::now();
        if (newCurve.mode == BERNSTEIN)
            newCurve.curvePoints = bezierCurve(newCurve.head, newCurve.steps);
        else
            newCurve.curvePoints = bezierCurveDeCasteljau(newCurve.controlPoints, newCurve.steps);
        auto end = std::chrono::high_resolution_clock::now();
        newCurve.computeTime = std::chrono::duration<double, std::milli>(end - start).count();
    
        if (bezierHead == NULL) {
            bezierHead = new ListChaineCourbes;
            bezierHead->list = newCurve;
            bezierHead->next = NULL;
        } else {
            ListChaineCourbes* current = bezierHead;
            while (current->next != NULL) {
                current = current->next;
            }
            current->next = new ListChaineCourbes;
            current->next->list = newCurve;
            current->next->next = NULL;
        }

        // Reset pour une nouvelle courbe
        currentControlPoints.clear();
        currentBezierMode = NONE;
    
        glutPostRedisplay();
    });
    glutAddMenuEntry("Tracer par Bernstein", 1);
    glutAddMenuEntry("Tracer par De Casteljau", 2);

    int BezierTestMenu = glutCreateMenu([](int value){
        if (value == 1) {
            currentControlPoints.clear();
            int N = 60; // par exemple 60 points de contrôle
            for (int i = 0; i < N; ++i) {
                Pixel p;
                p.x = 20 + (i * (windowWidth - 40)) / (N-1);
                p.y = 50 + (std::rand() % (windowHeight - 100)); // hauteur aléatoire
                p.color[0] = 1.0f; p.color[1] = 0.5f; p.color[2] = 0.0f;
                currentControlPoints.push_back(p);
            }
            drawingBezier = false;
            glutPostRedisplay();
            std::cout << "Test : " << N << " points de controle crees." << std::endl;
        }
    });
    glutAddMenuEntry("Générer 60 points Bézier (test perf)", 1);


    int menu = glutCreateMenu(menuCallBack);

    glutAddSubMenu("Couleurs", colorMenu);
    glutAddMenuEntry("polygone a decouper", 2);
    glutAddMenuEntry("trace fenetre", 3);
    glutAddSubMenu("fenetrage", FenMenu);
    glutAddSubMenu("remplissage", FillMenu);
    glutAddSubMenu("Poser points Bézier", PosePointsMenu);
    glutAddMenuEntry("Deplacer des points", 4);
    glutAddMenuEntry("Supprimer une courbe", 5);
    glutAddSubMenu("Tracer courbe Bézier", TraceBezierMenu);
    glutAddSubMenu("Test+50 point", BezierTestMenu);


    glutAttachMenu(GLUT_RIGHT_BUTTON);



    
    glutDisplayFunc(displayCallBack);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseCallBack);
    glutPassiveMotionFunc(passiveMotionCallback);
    glutMotionFunc(motionCallBack);
    glutReshapeFunc(reshapeCallback);

    glutTimerFunc(16, timerCallback, 0);

    glutMainLoop();

    return 0;
}