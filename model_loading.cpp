#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <learnopengl/filesystem.h>
#include <learnopengl/shader_m.h>
#include <learnopengl/camera.h>
#include <learnopengl/model.h>

#include <iostream>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 800;

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// camera
float radius = 5.0f;          // distance from model
float orbitYaw = 0.0f;        // horizontal angle (degrees)
float orbitPitch = 20.0f;     // vertical angle (degrees)
glm::vec3 target(20.0f, 200.0f, 0.0f);  // model center
glm::vec3 rewardposition(0.0f, 0.0f, 0.0f);  // model center

// --- Multiple Rocks ---
const int ROCK_COUNT = 3;
glm::vec3 rockPositions[ROCK_COUNT] = {
    glm::vec3(0.0f, 200.0f, 0.0f),
    glm::vec3(10.0f, 250.0f, 15.0f),
    glm::vec3(-15.0f, 220.0f, -5.0f)
};

glm::vec3 rockVelocities[ROCK_COUNT] = {
    glm::vec3(0.0f),
    glm::vec3(0.0f),
    glm::vec3(0.0f)
};


float targetYaw = orbitYaw;
float targetPitch = orbitPitch;

float smoothSpeed = 8.0f;  // higher = faster interpolation

glm::vec3 velocity(0.0f);   // vertical velocity
float gravity = -14.0f;     // gravity acceleration (units per second^2)
bool onGround = false;      // whether Wooloo is on the ground
bool isReleaseSpace = true;

struct TerrainBuffers {
    GLuint vao, vbo, ebo;
    GLsizei indexCount;
};

struct SphereCollider {
    glm::vec3 center;
    float radius;
};

struct TriangleCollider {
    glm::vec3 a, b, c;
};

void diamondSquare
(std::vector<std::vector<float>>& map, int size, float roughness)
{
    int stepSize = size - 1;
    float scale = roughness;

    while (stepSize > 1)
    {
        int halfStep = stepSize / 2;

        // ---- Diamond Step ----
        for (int x = 0; x < size - 1; x += stepSize)
        {
            for (int y = 0; y < size - 1; y += stepSize)
            {
                float avg = (map[x][y] +
                    map[x + stepSize][y] +
                    map[x][y + stepSize] +
                    map[x + stepSize][y + stepSize]) * 0.25f;

                map[x + halfStep][y + halfStep] = avg + ((rand() / (float)RAND_MAX) * 2 - 1) * scale;
            }
        }

        // ---- Square Step ----
        for (int x = 0; x < size; x += halfStep)
        {
            for (int y = (x + halfStep) % stepSize; y < size; y += stepSize)
            {
                float sum = 0.0f;
                int count = 0;

                if (x - halfStep >= 0) {
                    sum += map[x - halfStep][y];
                    count++;
                }
                if (x + halfStep < size) {
                    sum += map[x + halfStep][y];
                    count++;
                }
                if (y - halfStep >= 0) {
                    sum += map[x][y - halfStep];
                    count++;
                }
                if (y + halfStep < size) {
                    sum += map[x][y + halfStep];
                    count++;
                }

                map[x][y] = sum / count + ((rand() / (float)RAND_MAX) * 2 - 1) * scale;
            }
        }

        stepSize /= 2;
        scale *= 0.5f; // reduce displacement
    }
}

Mesh createTerrainMesh
(const std::vector<std::vector<float>>& map) {
    int size = map.size();
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    std::vector<glm::vec3> positions(size * size);
    std::vector<glm::vec3> normals(size * size, glm::vec3(0.0f));

    // --- Compute positions ---
    for (int z = 0; z < size; ++z) {
        for (int x = 0; x < size; ++x) {
            positions[z * size + x] = glm::vec3(
                (float)(x - size / 2), // X
                map[x][z],             // Y
                (float)(z - size / 2)  // Z
            );
        }
    }

    // --- Compute normals ---
    for (int z = 0; z < size - 1; ++z) {
        for (int x = 0; x < size - 1; ++x) {
            glm::vec3 p0 = positions[z * size + x];
            glm::vec3 p1 = positions[z * size + x + 1];
            glm::vec3 p2 = positions[(z + 1) * size + x];
            glm::vec3 n = glm::normalize(glm::cross(p1 - p0, p2 - p0));

            normals[z * size + x] += n;
            normals[z * size + x + 1] += n;
            normals[(z + 1) * size + x] += n;
        }
    }
    for (auto& n : normals) n = glm::normalize(n);

    // --- Create Vertex array ---
    for (int i = 0; i < size * size; ++i) {
        Vertex v;
        v.Position = positions[i];
        v.Normal = normals[i];
        v.TexCoords = glm::vec2((float)(i % size) / (size - 1), (float)(i / size) / (size - 1));
        v.Tangent = glm::vec3(0.0f);
        v.Bitangent = glm::vec3(0.0f);
        std::fill(std::begin(v.m_BoneIDs), std::end(v.m_BoneIDs), 0);
        std::fill(std::begin(v.m_Weights), std::end(v.m_Weights), 0.0f);
        vertices.push_back(v);
    }

    // --- Create indices for triangle strip ---
    for (int z = 0; z < size - 1; ++z) {
        for (int x = 0; x < size - 1; ++x) {
            int topLeft = z * size + x;
            int topRight = topLeft + 1;
            int bottomLeft = (z + 1) * size + x;
            int bottomRight = bottomLeft + 1;

            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    std::vector<Texture> textures;

    // Suppose you have a file path to a texture
    std::string texturePath = FileSystem::getPath("resources/objects/Pokemon/RockTexture3.png");

    //std::cout << texturePath << std::endl;
    // Load the texture using your existing function
    unsigned int textureID = TextureFromFile(texturePath.c_str(), "", false);

    // Create a Texture struct and fill its fields
    Texture grassTexture;
    grassTexture.id = textureID;
    grassTexture.type = "texture_diffuse"; // match your shader naming convention
    grassTexture.path = texturePath;

    // Add to the vector
    textures.push_back(grassTexture);

    return Mesh(vertices, indices, textures);
}

float getTerrainHeight
(const std::vector<std::vector<float>>& heightmap, float x, float z)
{
    int size = heightmap.size();
    int half = size / 2;

    // Convert world (x, z) to heightmap grid coordinates
    float fx = x + half;
    float fz = z + half;

    int ix = (int)floor(fx);
    int iz = (int)floor(fz);

    if (ix < 0 || iz < 0 || ix >= size - 1 || iz >= size - 1)
        return 0.0f; // out of bounds

    // Bilinear interpolation
    float tx = fx - ix;
    float tz = fz - iz;

    float h00 = heightmap[ix][iz];
    float h10 = heightmap[ix + 1][iz];
    float h01 = heightmap[ix][iz + 1];
    float h11 = heightmap[ix + 1][iz + 1];

    float h0 = glm::mix(h00, h10, tx);
    float h1 = glm::mix(h01, h11, tx);

    return glm::mix(h0, h1, tz);
}

void renderTerrain
(GLuint vbo, int size)
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexPointer(3, GL_FLOAT, 0, 0);

    glColor3f(0.3f, 1.0f, 0.3f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    for (int z = 0; z < size - 1; ++z)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int x = 0; x < size; ++x)
        {
            glArrayElement((z + 1) * size + x);
            glArrayElement(z * size + x);
        }
        glEnd();
    }

    glDisableClientState(GL_VERTEX_ARRAY);
}

glm::vec3 ClosestPtPointTriangle
(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c)
{
    glm::vec3 ab = b - a;
    glm::vec3 ac = c - a;
    glm::vec3 ap = p - a;

    float d1 = glm::dot(ab, ap);
    float d2 = glm::dot(ac, ap);

    if (d1 <= 0.0f && d2 <= 0.0f) return a; // barycentric (1,0,0)

    glm::vec3 bp = p - b;
    float d3 = glm::dot(ab, bp);
    float d4 = glm::dot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) return b; // barycentric (0,1,0)

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        float v = d1 / (d1 - d3);
        return a + v * ab; // barycentric (1−v, v, 0)
    }

    glm::vec3 cp = p - c;
    float d5 = glm::dot(ab, cp);
    float d6 = glm::dot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) return c; // barycentric (0,0,1)

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        float w = d2 / (d2 - d6);
        return a + w * ac; // barycentric (1−w, 0, w)
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric (0, 1−w, w)
    }

    // inside face region
    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

bool TestSphereTriangle
(const SphereCollider& sphere, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, glm::vec3& closestPoint)
{
    closestPoint = ClosestPtPointTriangle(sphere.center, a, b, c);
    glm::vec3 v = closestPoint - sphere.center;
    return glm::dot(v, v) <= sphere.radius * sphere.radius;
}

int main()
{
    srand(time(NULL));
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    
    glfwSetWindowPos(window, 400, 200);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // tell stb_image.h to flip loaded texture's on the y-axis (before loading model).
    stbi_set_flip_vertically_on_load(false);

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile shaders
    // -------------------------
    Shader ourShader("1.model_loading.vs", "1.model_loading.fs");

    // load models
    // -----------
    Model charizard(FileSystem::getPath("resources/objects/Pokemon/charizard.obj"));
    Model wooloo(FileSystem::getPath("resources/objects/Pokemon/wooloo.obj"));
    Model pokeball(FileSystem::getPath("resources/objects/Pokemon/pokeball.obj"));
    Model rock(FileSystem::getPath("resources/objects/Pokemon/rock.obj"));

    // Generate heightmap
    int n = 6; // produces 2^7 + 1 = 129 x 129 grid
    int size = (1 << n) + 1;
    
    std::vector<std::vector<float>> heightmap(size, std::vector<float>(size, 0.0f));
    diamondSquare(heightmap, size, 10.0f);

    Mesh terrainMesh = createTerrainMesh(heightmap);

    rewardposition = glm::vec3(0.0f, getTerrainHeight(heightmap, 0.0f, 0.0f) + 1.0f, 0.0f);
    // draw in wireframe
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        //std::cout << "velocity.y = " << velocity.y << std::endl;
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // don't forget to enable shader before setting uniforms
        ourShader.use();

        //// view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        //glm::mat4 view = camera.GetViewMatrix();
        ourShader.setMat4("projection", projection);
        //ourShader.setMat4("view", view);
        
        // Smoothly interpolate camera orbit
        float lerpFactor = 1.0f - expf(-smoothSpeed * deltaTime);
        orbitYaw = glm::mix(orbitYaw, targetYaw, lerpFactor);
        orbitPitch = glm::mix(orbitPitch, targetPitch, lerpFactor);

        // Camera uses smoothed orbit
        float smoothYawRad = glm::radians(orbitYaw);
        float smoothPitchRad = glm::radians(orbitPitch);

        // Model uses target orbit (instant)
        float targetYawRad = glm::radians(targetYaw);
        float targetPitchRad = glm::radians(targetPitch);

        // Calculate camera position (smooth)
        float camX = radius * cos(smoothPitchRad) * sin(smoothYawRad);
        float camY = radius * sin(smoothPitchRad);
        float camZ = radius * cos(smoothPitchRad) * cos(smoothYawRad);

        glm::vec3 cameraPos = target + glm::vec3(camX, camY, camZ);
        glm::mat4 view = glm::lookAt(cameraPos, target, glm::vec3(0.0f, 1.0f, 0.0f));
        ourShader.setMat4("view", view);

        // Make wooloo face the same direction as the camera
        glm::vec3 modelDir;

        targetPitchRad = 0;
        modelDir.x = cos(targetPitchRad) * sin(targetYawRad);
        modelDir.y = sin(targetPitchRad);
        modelDir.z = cos(targetPitchRad) * cos(targetYawRad);
        glm::vec3 cameraDir = -glm::normalize(modelDir);
        glm::vec3 worldUp(0.0f, 1.0f, 0.0f);

        // Construct a rotation matrix that faces the same direction
        glm::vec3 right = glm::normalize(glm::cross(worldUp, cameraDir));
        glm::vec3 up = glm::normalize(glm::cross(cameraDir, right));

        glm::mat4 rotation = glm::mat4(1.0f);
        rotation[0] = glm::vec4(right, 0.0f);
        rotation[1] = glm::vec4(up, 0.0f);
        rotation[2] = glm::vec4(cameraDir, 0.0f);

        // render the loaded model
        glm::mat4 woolooModel = glm::mat4(1.0f);
        glm::vec3 woolooPosition = glm::vec3(target.x, target.y - 10.0f, target.z);
        woolooModel = glm::translate(woolooModel, target);
        woolooModel *= rotation; // apply facing rotation
        woolooModel = glm::scale(woolooModel, glm::vec3(1.0f));
        ourShader.setMat4("model", woolooModel);
        wooloo.Draw(ourShader);


        for (int i = 0; i < ROCK_COUNT; ++i)
        {
            glm::mat4 rockModel = glm::mat4(1.0f);
            rockModel = glm::translate(rockModel, rockPositions[i]);
            rockModel = glm::scale(rockModel, glm::vec3(3.0f));
            ourShader.setMat4("model", rockModel);
            rock.Draw(ourShader);
        }

        glm::mat4 pokeballModel = glm::mat4(1.0f);
        pokeballModel = glm::translate(pokeballModel, rewardposition); // translate it down so it's at the center of the scene
        pokeballModel = glm::scale(pokeballModel, glm::vec3(0.1f, 0.1f, 0.1f));	// it's a bit too big for our scene, so scale it down
        ourShader.setMat4("model", pokeballModel);
        pokeball.Draw(ourShader);

        glm::mat4 terrainModel = glm::mat4(1.0f);
        ourShader.setMat4("model", terrainModel);
        terrainMesh.Draw(ourShader);


        SphereCollider woolooCollider;
        woolooCollider.radius = 0.2f; // adjust based on model size

        SphereCollider pokeballCollider;
        pokeballCollider.radius = 2.0f; // adjust based on model size

        SphereCollider rockCollider;
        rockCollider.radius = 3.0f; // adjust based on model size

        // --- Gravity ---
        if (!onGround) {
            velocity.y += gravity * deltaTime;  // v = v + g*dt
           
        }
        else if (onGround && velocity.y < 0.0f) {
            velocity.y = 0.0f; // stop vertical motion if on ground
            
        }

        // Apply vertical velocity
        target.y += velocity.y * deltaTime;

        for (int i = 0; i < ROCK_COUNT; ++i)
        {
            rockVelocities[i] += glm::vec3(0.0f, gravity * deltaTime, 0.0f);
            rockPositions[i] += rockVelocities[i] * deltaTime;

            // Reset rock if it falls off
            if (rockPositions[i].y <= -20.0f)
            {
                rockPositions[i] = glm::vec3(
                    (rand() % 100 - 50) * 0.5f,
                    200.0f + rand() % 50,
                    (rand() % 100 - 50) * 0.5f
                );
                rockVelocities[i] = glm::vec3(0.0f);
            }
        }

        // Check wooloo-pokeball colliding
        woolooCollider.center = target;
        // Update pokeball collider
        pokeballCollider.center = rewardposition;

        // --- Wooloo and Pokeball Collision ---
        glm::vec3 diff = woolooCollider.center - pokeballCollider.center;
        float dist2 = glm::dot(diff, diff);
        float rSum = woolooCollider.radius + pokeballCollider.radius;

        if (dist2 <= rSum * rSum)
        {
            std::cout << "Wooloo collected Pokeball!" << std::endl;

            float newX = (rand() % 100 - 50) * 0.5f;
            float newZ = (rand() % 100 - 50) * 0.5f;
            float newY = getTerrainHeight(heightmap, newX, newZ) + 2.0f;

            rewardposition = glm::vec3(newX, newY, newZ);
        }

        // Check wooloo-terrain colliding
        onGround = false;

        for (int i = 0; i < terrainMesh.indices.size(); i += 3)
        {
            glm::vec3 a = terrainMesh.vertices[terrainMesh.indices[i]].Position;
            glm::vec3 b = terrainMesh.vertices[terrainMesh.indices[i + 1]].Position;
            glm::vec3 c = terrainMesh.vertices[terrainMesh.indices[i + 2]].Position;

            glm::vec3 closest;
            if (TestSphereTriangle(woolooCollider, a, b, c, closest))
            {
                // Collision response: push Wooloo above terrain
                glm::vec3 n = glm::normalize(glm::cross(b - a, c - a));
                float penetration = woolooCollider.radius - glm::length(closest - woolooCollider.center);
                target += n * penetration;
                woolooCollider.center = target;

                // Check if normal points roughly up => consider on ground
                if (n.y > 0.5f) onGround = true;

                break; // stop after first collision
            }

            
        }

        // Check rock-terrain colliding
        float bounceFactor = 0.4f;

        for (int r = 0; r < ROCK_COUNT; ++r)
        {
            SphereCollider rockCollider;
            rockCollider.center = rockPositions[r];
            rockCollider.radius = 2.0f;

            for (int i = 0; i < terrainMesh.indices.size(); i += 3)
            {
                glm::vec3 a = terrainMesh.vertices[terrainMesh.indices[i]].Position;
                glm::vec3 b = terrainMesh.vertices[terrainMesh.indices[i + 1]].Position;
                glm::vec3 c = terrainMesh.vertices[terrainMesh.indices[i + 2]].Position;

                glm::vec3 closest;
                if (TestSphereTriangle(rockCollider, a, b, c, closest))
                {
                    glm::vec3 n = glm::normalize(glm::cross(b - a, c - a));
                    float penetration = rockCollider.radius - glm::length(closest - rockCollider.center);
                    rockPositions[r] += n * penetration;

                    // Bounce
                    rockVelocities[r] = glm::reflect(rockVelocities[r], n) * bounceFactor;
                    if (rockVelocities[r].y < 0.2f) { rockVelocities[r].y = 12.0f; }

                    // Small velocity cutoff
                    if (glm::length(rockVelocities[r]) < 0.1f)
                        rockVelocities[r] = glm::vec3(0.0f);

                    break;
                }
            }
        }

        // CHeck wooloo-rock colliding
        for (int r = 0; r < ROCK_COUNT; ++r)
        {
            SphereCollider rockCollider;
            rockCollider.center = rockPositions[r];
            rockCollider.radius = 2.0f;

            glm::vec3 diff = woolooCollider.center - rockCollider.center;
            float dist = glm::length(diff);
            float minDist = woolooCollider.radius + rockCollider.radius;

            if (dist < minDist && dist > 0.0001f)
            {
                glm::vec3 pushDir = glm::normalize(diff);
                float penetration = minDist - dist;

                // Push Wooloo away
                target += pushDir * penetration * 0.5f;

                // Apply a small impulse to Wooloo
                velocity += pushDir * 5.0f;

                // Also nudge rock slightly backward
                rockPositions[r] -= pushDir * penetration * 0.5f;
            }
        }


        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);


    float moveSpeed = 5.0f;

    glm::vec3 camForward = glm::normalize(glm::vec3(
        -cos(glm::radians(targetPitch)) * sin(glm::radians(targetYaw)),
        0.0f,
        -cos(glm::radians(targetPitch)) * cos(glm::radians(targetYaw))
    ));
    glm::vec3 camRight = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), camForward));


    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        target += camForward * moveSpeed * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        target -= camForward * moveSpeed * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        target += camRight * moveSpeed * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        target -= camRight * moveSpeed * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
    {
        if (isReleaseSpace == true && onGround == true)
        {
            //std::cout << "Press space" << std::endl;

            velocity.y = 10.0f;

            isReleaseSpace = false;
        }
        
    }
    else
    {
        //std::cout << " " << std::endl;
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_RELEASE)
    {
        isReleaseSpace = true;
    }

        
        
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{

    static float lastX = SCR_WIDTH / 2.0f;
    static float lastY = SCR_HEIGHT / 2.0f;
    static bool firstMouse = true;

    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    // Apply to target rotation only
    targetYaw -= xoffset;
    targetPitch -= yoffset;

    // Clamp pitch
    if (targetPitch > 89.0f) targetPitch = 89.0f;
    if (targetPitch < -89.0f) targetPitch = -89.0f;
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    radius -= static_cast<float>(yoffset);
    if (radius < 2.0f) radius = 2.0f;
    if (radius > 20.0f) radius = 20.0f;

    //camera.ProcessMouseScroll(static_cast<float>(yoffset));
}
