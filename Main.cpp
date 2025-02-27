#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <stack>
#include <cmath>


enum BorderType { LEFT, RIGHT, BOTTOM, TOP };

struct Polygons {
    std::vector<std::pair<int, int>> points;
    float color[3];
};

struct SelectionBox {
    std::pair<int, int> start;
    std::pair<int, int> end;    
    float color[3];
};

float selectedColor[3] = {1.0f, 1.0f, 1.0f};
int mouseX = 0, mouseY = 0;
int windowWidth = 500, windowHeight = 500;
bool starSelectionBox = false;
bool showPixel = false;
bool drawTrail = false;
bool showLine = false;
bool drawingPolygon = false;
bool drawingSelectionBox = false;
bool makeFill = false;
std::vector<SelectionBox> selectionBoxList;
std::pair<int, int> selectionBoxStart;
std::pair<int, int> selectionBoxEnd;
std::vector<std::pair<int, int>> trail;
std::vector<std::pair<int, int>> LinePoints;
std::vector<Polygons> polygonList;
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
    // On lit 1 pixel en format RGB8
    glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, px);
    // Convertit [0..255] => [0..1]
    outColor[0] = px[0] / 255.0f;
    outColor[1] = px[1] / 255.0f;
    outColor[2] = px[2] / 255.0f;
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
    // Pour éviter d'attendre la fin du dessin, on peut forcer le flush
    // glFlush(); // si vous êtes en simple buffer
    // ou glutSwapBuffers(); // si vous êtes en double buffer
    // Selon votre code, vous pouvez laisser glutSwapBuffers()
    // se faire dans displayCallBack().
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


void floodFillRecursive(const std::pair<int,int>& P,float boundaryColor[3] ,float newColor[3])
{
    int x = P.first;
    int y = P.second;

    // 1. Vérification des coordonnées
    if (x < 0 || x >= windowWidth || y < 0 || y >= windowHeight) {
        return; 
    }

    // 2. Lire la couleur actuelle du pixel
    float currentColor[3];
    getPixelColor(x, y, currentColor);

    // 3. Condition d'arrêt:
    //    - si c'est la frontière
    //    - ou si c'est déjà la couleur de remplissage
    if (isSameColor(currentColor, boundaryColor) ||
        isSameColor(currentColor, newColor)) {
        return; 
    }

    // 4. Colorier le pixel
    putPixel(x, y, newColor);

    // 5. Appels récursifs dans les 4 directions
    floodFillRecursive(std::make_pair(x+1, y), boundaryColor, newColor);
    floodFillRecursive(std::make_pair(x-1, y), boundaryColor, newColor);
    floodFillRecursive(std::make_pair(x, y+1), boundaryColor, newColor);
    floodFillRecursive(std::make_pair(x, y-1), boundaryColor, newColor);

    // Si vous désirez un 8-voisinage:
    // floodFillRecursive(std::make_pair(x+1, y+1), newColor);
    // floodFillRecursive(std::make_pair(x-1, y-1), newColor);
    // floodFillRecursive(std::make_pair(x+1, y-1), newColor);
    // floodFillRecursive(std::make_pair(x-1, y+1), newColor);
}


std::pair<int,int> computeIntersection(
    const std::pair<int,int>& S1, 
    const std::pair<int,int>& S2, 
    float boundaryValue, 
    BorderType border
) {
    float x1 = (float)S1.first,  y1 = (float)S1.second;
    float x2 = (float)S2.first,  y2 = (float)S2.second;
    float x, y;
    
    if (border == LEFT || border == RIGHT) {
        // Eviter la division par zéro
        if (std::abs(x2 - x1) < 1e-5f) {
            // Segment quasi vertical -> On peut retourner S1 (ou S2) 
            // ou gérer autrement
            return S1; 
        }
        float t = (boundaryValue - x1) / (x2 - x1);
        x = boundaryValue;
        y = y1 + t * (y2 - y1);
    } 
    else {
        // border == TOP ou BOTTOM (clip horizontal)
        if (std::abs(y2 - y1) < 1e-5f) {
            return S1;
        }
        float t = (boundaryValue - y1) / (y2 - y1);
        y = boundaryValue;
        x = x1 + t * (x2 - x1);
    }

    return std::make_pair((int)x, (int)y);
}

std::vector<std::pair<int,int>> clipWithBorder(
    const std::vector<std::pair<int,int>>& inPolygon,
    float boundaryValue,
    BorderType border
) {
    std::vector<std::pair<int,int>> outPolygon;
    // On "ferme" le polygone : dernier -> premier
    for (size_t i = 0; i < inPolygon.size(); i++) {
        std::pair<int,int> current = inPolygon[i];
        std::pair<int,int> next = inPolygon[(i + 1) % inPolygon.size()];

        bool currInside = isInside(current, boundaryValue, border);
        bool nextInside = isInside(next, boundaryValue, border);

        if (currInside && nextInside) {
            // Garder le next
            outPolygon.push_back(next);
        } 
        else if (currInside && !nextInside) {
            // De dedans à dehors => intersection
            outPolygon.push_back(computeIntersection(current, next, boundaryValue, border));
        }
        else if (!currInside && nextInside) {
            // De dehors à dedans => intersection + next
            outPolygon.push_back(computeIntersection(current, next, boundaryValue, border));
            outPolygon.push_back(next);
        }
        // Sinon dehors->dehors => rien
    }
    return outPolygon;
}

Polygons clipPolygon(const Polygons& inputPoly, const SelectionBox& box) {
    Polygons output;

    // Déterminer xmin, xmax, ymin, ymax
    float xmin = std::min(box.start.first, box.end.first);
    float xmax = std::max(box.start.first, box.end.first);
    float ymin = std::min(box.start.second, box.end.second);
    float ymax = std::max(box.start.second, box.end.second);

    // On récupère les points d'entrée
    std::vector<std::pair<int,int>> clippedPoints = inputPoly.points;

    // 1) Clip bord gauche : x >= xmin
    clippedPoints = clipWithBorder(clippedPoints, xmin, LEFT);

    // 2) Clip bord droit : x <= xmax
    clippedPoints = clipWithBorder(clippedPoints, xmax, RIGHT);

    // 3) Clip bord bas : y >= ymin
    clippedPoints = clipWithBorder(clippedPoints, ymin, BOTTOM);

    // 4) Clip bord haut : y <= ymax
    clippedPoints = clipWithBorder(clippedPoints, ymax, TOP);

    // Assigner le nouveau polygone
    output.points = clippedPoints;
    
    // Vous pouvez choisir une couleur différente
    output.color[0] = 1.0f; // Magenta, par ex.
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


void menuCallBack(int value){
    std::vector<Polygons> newPolygon;
    bool clippedFlag = false;
    switch (value) {
        case 1:
            break;
        case 2:
            std::cout << "Option 2 selectionnee" << std::endl;
            showPixel = true;
            showLine = true;
            break;
        case 3:
            std::cout << "Option 3 selectionnee" << std::endl;
            drawingSelectionBox = true;
            break;
        case 4:
            std::cout << "Option 4 selectionnee" << std::endl;
            if (!polygonList.empty() && !selectionBoxList.empty()) {
                for(Polygons p : polygonList){
                    for(SelectionBox sb : selectionBoxList){

                        Polygons clipped = clipPolygon(p, sb);
                        // On l'ajoute à la liste
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
        case 5:
            std::cout << "Option 5 selectionnee" << std::endl;
            makeFill = true;
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
         // Ajouter la position actuelle à la traînée
        if (showLine && showPixel) {
                drawTrail = true;
                trail.emplace_back(x, windowHeight - y);
                LinePoints.emplace_back(x, windowHeight - y); // Ajouter un point au polygone
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
        if (makeFill && !drawingPolygon) {
        int adjustedY = windowHeight - y; // pour ajuster la coordonnée y
            // Parcourir les polygones terminés
            for (auto& poly : polygonList) {
                if (pointInPolygon(x, adjustedY, poly.points)) {
                    std::cout << "Remplissage du polygone à partir du point (" 
                            << x << ", " << adjustedY << ")" << std::endl;
                    // On peut utiliser ici la version récursive ou itérative
                    floodFillRecursive(std::make_pair(x, adjustedY), poly.color, selectedColor);
                    // Une fois le remplissage effectué, on réinitialise le flag
                    makeFill = false;
                    break;
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
            polygonList.push_back(currentPolygon); // Ajouter le polygone actuel à la liste
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
            break;
    }
    glutPostRedisplay();
}

void timerCallback(int value) {
    glutPostRedisplay(); // Forcer un redessin
    glutTimerFunc(16, timerCallback, 0); // Réenregistrer le timer pour 16 ms (environ 60 FPS)
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

    int menu = glutCreateMenu(menuCallBack);

    glutAddSubMenu("Couleurs", colorMenu);
    glutAddMenuEntry("polygone a decouper", 2);
    glutAddMenuEntry("trace fenetre", 3);
    glutAddMenuEntry("fenetrage", 4);
    glutAddMenuEntry("remplissage", 5);


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