#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "main.h"
#include "dosusangWindow.c"
using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth = 255;
float* zbuffer = new float[width * height];

void line(Vec2i t0, Vec2i t1, TGAImage& image, TGAColor color) {
    int x0 = t0.x;
    int x1 = t1.x;

    int y0 = t0.y;
    int y1 = t1.y;

    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2;
    int error2 = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) {
            y += (y1 > y0 ? 1 : -1);
            error2 -= dx * 2;
        }
    }
}

void drawFrame(TGAImage &image) {
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        for (int j = 0; j < 3; j++) {
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j + 1) % 3]);
            int x0 = (v0.x + 1.) * width / 2.;
            int y0 = (v0.y + 1.) * height / 2.;
            int x1 = (v1.x + 1.) * width / 2.;
            int y1 = (v1.y + 1.) * height / 2.;
            line(Vec2i(x0, y0),Vec2i(x1, y1), image, white);
        }
    }
}


// 返回p点对于三角形的质心坐标
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f AB = B - A;
    Vec3f AC = C - A;
    Vec3f PA = A - P;
    Vec3f tar = Vec3f(AB.x, AC.x, PA.x) ^ Vec3f(AB.y, AC.y, PA.y);

    if (abs(tar.z) < 1e-2) return Vec3f(-1, 1, 1);
    return Vec3f(1.0f - (tar.x + tar.y) / tar.z, tar.x / tar.z, tar.y / tar.z);
}

void triangle(Vec3f t0, Vec3f t1, Vec3f t2, TGAImage& img, TGAColor color) {
    int minX = max(min(min(t0.x, t1.x), t2.x), 0.0f);
    int minY = max(min(min(t0.y, t1.y), t2.y), 0.0f);
    int maxX = min(max(max(t0.x, t1.x), t2.x), width-1.0f);
    int maxY = min(max(max(t0.y, t1.y), t2.y), height-1.0f);
    Vec3f P;
    for (int i = minX; i <= maxX; i++) {
        for (int j = minY; j <= maxY; j++) {
            P.x = i;
            P.y = j;
            Vec3f bc_screen = barycentric(t0, t1, t2, P);
            // 三角形剔除
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            // 深度剔除
            P.z = 0;
            P.z = bc_screen * Vec3f(t0.z, t1.z, t2.z);
            if (zbuffer[int(P.x + P.y * width)] < P.z) {
               zbuffer[int(P.x + P.y * width)] = P.z;
                img.set(P.x, P.y, color);
            }

        }
    }
}

void triangle(Vec3f t0, Vec3f t1, Vec3f t2, TGAImage& img, Vec2i uv0, Vec2i uv1, Vec2i uv2) {
    int minX = max(min(min(t0.x, t1.x), t2.x), 0.0f);
    int minY = max(min(min(t0.y, t1.y), t2.y), 0.0f);
    int maxX = min(max(max(t0.x, t1.x), t2.x), width - 1.0f);
    int maxY = min(max(max(t0.y, t1.y), t2.y), height - 1.0f);
    Vec3f P;
    for (int i = minX; i <= maxX; i++) {
        for (int j = minY; j <= maxY; j++) {
            P.x = i;
            P.y = j;
            Vec3f bc_screen = barycentric(t0, t1, t2, P);
            // 三角形剔除
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            // 深度剔除
            P.z = 0;
            // 深度插值
            P.z = bc_screen * Vec3f(t0.z, t1.z, t2.z);
            // UV插值
            Vec2i uv = uv0 * bc_screen.x + uv1 * bc_screen.y + uv2 * bc_screen.z;
            if (zbuffer[int(P.x + P.y * width)] < P.z) {
                zbuffer[int(P.x + P.y * width)] = P.z;
                img.set(P.x, P.y, model->diffuse(uv));
            }
        }
    }
}

// 视口变换
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x + w / 2.f;
    m[1][3] = y + h / 2.f;
    m[2][3] = depth / 2.f;

    m[0][0] = w / 2.f;
    m[1][1] = h / 2.f;
    m[2][2] = depth / 2.f;
    return m;
}

// View矩阵，转换到相机空间
Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye - center).normalize();
    Vec3f x = (up ^ z).normalize();
    Vec3f y = (z ^ x).normalize();
    Matrix res = Matrix::identity(4);
    for (int i = 0; i < 3; i++) {
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    TGAImage image(width, height, TGAImage::RGB);

    // 初始化Zbuffer
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    // 定义一个方向光
    Vec3f light_dir(0, 0, -1); 

    Vec3f eye(1, 10, 3);
    Vec3f center(0, 0, 0);

    // 世界坐标转到相机坐标的矩阵 平移+旋转
    Matrix ModelView = lookat(eye, center, Vec3f(0, 1, 0));

    // 投影矩阵 投影
    Matrix Projection = Matrix::identity(4);
    // 这里投影平面是相机空间下的模型原点所处的xoy平面，缩放因子就是-1/dist(eye, center)
    Projection[3][2] = -1.f / (eye - center).norm();

    // 视口变换 缩放
    Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        Vec2i model_uv[3];

        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec3f(ViewPort * Projection * ModelView * Matrix(v));
            world_coords[j] = v;
            model_uv[j] = model->uv(i, j);
        }
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;

        triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, model_uv[0], model_uv[1], model_uv[2]);
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;

    return 0;
}

