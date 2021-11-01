#include "Operators.cpp"
#include "HalfedgeBuilder.cpp"


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
using namespace Eigen;




bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == '1')
    {
        viewer.data().clear();
        
        
    }
    else if (key == '2')
    {
        viewer.data().clear();
        
    }
    return false;
}

MatrixXd areaMatrix(MatrixXd& V, MatrixXi& F) {
    MatrixXd res = MatrixXd::Zero(F.rows(), 1);
    Vector3d v1(0.0,0.0,0.0);
    Vector3d v2(0.0,0.0,0.0);
    Vector3d v3(0.0,0.0,0.0);
    Vector3d a, b;

    for (int i = 0; i < F.rows(); ++i) {
        v1 << V(F(i, 0), 0), V(F(i, 0), 1), V(F(i, 0), 2);
        v2 << V(F(i, 1), 0), V(F(i, 1), 1), V(F(i, 1), 2);
        v3 << V(F(i, 2), 0), V(F(i, 2), 1), V(F(i, 2), 2);
        a = v2 - v1;
        b = v3 - v1;
        res(i, 0) = (1.0 / 2) * (a.cross(b)).norm();
    }
    return res;
}
//this method will calculate the neighbour faces for a given vertex v
vector<int> vertexNeighbour(HalfedgeDS he, int v) {
    vector <int> n;
    int edge = he.getEdge(v);
    int oppEdge = he.getOpposite(edge);
    int prevEdge = he.getPrev(oppEdge);
    if (he.getFace(prevEdge) != -1) {
        n.push_back(he.getFace(prevEdge));
    }
    while (prevEdge != edge) {
        oppEdge = he.getOpposite(prevEdge);
        prevEdge = he.getPrev(oppEdge);
        if (he.getFace(prevEdge) != -1) {
            n.push_back(he.getFace(prevEdge));
        }
    }
    return n;
}



bool elementexists(vector<int> v, int j) {
    for (int i = 0; i < v.size(); i = i + 1) {
        if (v.at(i) == j) return true;
    }
    return false;
}

//this method will calculate the neighbour faces for a given face
vector<int> faceNeighbours(HalfedgeDS he, MatrixXd &V, MatrixXi& F,int f) {
    vector<int> res;
    for (int j = 0; j < 3; j = j + 1) {
            vector<int> a = vertexNeighbour(he, F(f, j));
            for (int i = 0; i < a.size(); i = i + 1) {
                if (!(elementexists(res, a.at(i)))) {
                    res.push_back(a.at(i));
                }
            }
        }
    return res;
}

int nextTriangleInColorFlow(HalfedgeDS he, MatrixXd& V, MatrixXi& F, int f) {
    return 0;
}


void buildNeighbours(std::map<int, vector<int>>& adj, HalfedgeDS he, MatrixXd& V) {
    for (int i = 0; i < V.rows(); i++) {
        adj.insert({ i, vertexNeighbour(he, i) });
    }
}




MatrixXd divMatrix1(std::map<int, vector<int>>& adj, MatrixXd& V,MatrixXi &F, MatrixXd& areaMatrix,MatrixXd & massMatrix, MatrixXd& grad) {
    MatrixXd div = MatrixXd::Zero(V.rows(), F.rows() * 3);

    for (int i = 0; i < V.rows(); i++) {
        
        for (const auto& value : adj.at(i)) {
                div(i, 3 * value) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value, i);
                div(i, 3 * value + 1) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value + 1, i);
                div(i, 3 * value + 2) = -(areaMatrix(value, 0) / massMatrix(i, i)) * grad(3 * value + 2, i);
        }
    }
    
    return div;
}

double distance(MatrixXd &V1, MatrixXd &V2) {
    return (V2 - V1).norm();
}

double cosinus_law(double a, double b, double c) {
    return acos((-pow(c, 2) + pow(a, 2) +pow(b, 2)) / (2.0 * a * b));
}

double cotan(double i) { return(1 / tan(i)); }

MatrixXd internal_angles(MatrixXd& V, MatrixXi& F) {
    MatrixXd angles(F.rows(), 3);
    double a=0, b=0, c=0;
    MatrixXd A(1,3), B(1,3), C(1,3);
    for (int i = 0; i < F.rows(); i = i + 1) {
        A = V.row(F(i, 0));
        B = V.row(F(i, 1));
        C = V.row(F(i, 2));

        a = distance(C, B);
        b = distance(A, C);
        c = distance(A, B);

        

        angles(i, 0) = cosinus_law(c, b, a);
        angles(i, 1) = cosinus_law(a, c, b);
        angles(i, 2) = cosinus_law(a, b, c);

    }

    return angles;
}

double cross(MatrixXd& e, MatrixXd& X) {
    Vector3d e1(e(0, 0), e(0, 1), e(0, 2));
    Vector3d X1(X(0, 0), X(0, 1), X(0, 2));
    return (e1.dot(X1));
}



MatrixXd divMatrix(std::map<int, vector<int>>& adj, MatrixXd& V, MatrixXi& F, MatrixXd& areaMatrix, MatrixXd& massMatrix, MatrixXd& gradU) {
    MatrixXd div = MatrixXd::Zero(V.rows(),1);
    MatrixXd angles = internal_angles(V, F);
   
    MatrixXd X(1, 3), e1(1, 3), e2(1, 3);
    vector<int> e,ind;
    
    for (int i = 0; i < V.rows(); i++) {
        
        for (const auto& value : adj.at(i)) {
            X = gradU.row(value);
            for (int j = 0; j < 3; j = j + 1) {
                if (F(value, j) != i) 
                    e.push_back(F(value, j));
                    ind.push_back(j);
            }
            e1 = V.row(e.at(0)) - V.row(i);
            e2 = V.row(e.at(1)) - V.row(i);
            div(i, 0) += 0.5*cotan(angles(value,ind.at(0)))*cross(e1,X)+0.5*cotan(angles(value,ind.at(1)))*cross(e2,X);
            
            e.clear();
            ind.clear();
        }
        }
    

    return div;
}
void set_colormap(igl::opengl::glfw::Viewer& viewer)
{
    const int num_intervals = 30;
    Eigen::MatrixXd CM(num_intervals, 3);
    // Colormap texture
    for (int i = 0; i < num_intervals; i++)
    {
        double t = double(num_intervals - i - 1) / double(num_intervals - 1);
        CM(i, 0) = std::max(std::min(2.0 * t - 0.0, 1.0), 0.0);
        CM(i, 1) = std::max(std::min(2.0 * t - 1.0, 1.0), 0.0);
        CM(i, 2) = std::max(std::min(6.0 * t - 5.0, 1.0), 0.0);
    }
    igl::isolines_map(Eigen::MatrixXd(CM), CM);
    viewer.data().set_colormap(CM);
}



void showVectorField(MatrixXd& V, MatrixXi& F, MatrixXd& GU, igl::opengl::glfw::Viewer& viewer) {
    // Draw a black segment in direction of gradient at face barycenters
    MatrixXd BC;
    igl::barycenter(V, F, BC);
    const RowVector3d red(1, 0, 0);
    const RowVector3d blue(0, 0, 1);
    viewer.data().add_edges(BC, BC + GU * 0.5, red);
    viewer.data().add_edges(BC + GU * 0.5, BC + GU * 0.5 + GU * 0.5, blue);
}

void showMaxMinCurvature(MatrixXd& V, MatrixXi& F, igl::opengl::glfw::Viewer& viewer) {
    MatrixXd PD1, PD2;
    VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    // mean curvature
    MatrixXd H = 0.5 * (PV1 + PV2);


    viewer.data().set_data(H);

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V, F);

    // Draw a red segment parallel to the maximal curvature direction
    const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
    viewer.data().add_edges(V + PD1 * avg, V - PD1 * avg, red);

    // Draw a blue segment parallel to the minimal curvature direction
    viewer.data().add_edges(V + PD2 * avg, V - PD2 * avg, blue);
}

int main()
{
    // Load a mesh in OFF format
    igl::opengl::glfw::Viewer viewer;
    Eigen::MatrixXd V, C;
    Eigen::MatrixXi F;

    Eigen::MatrixXd I,A;

    std::string filename = "C:\\Users\\Rached\\Documents\\Telecom\\igd\\test\\Mesh Coloring\\Mesh Coloring\\data\\sphere.off";
    igl::read_triangle_mesh(filename, V, F);
    
    HalfedgeBuilder* builder = new HalfedgeBuilder();  

    HalfedgeDS he = builder->createMesh(V.rows(), F);  
    
    std::map<int, vector<int>> adj;
    Vector3d qsdqs = V.row(0);

    buildNeighbours(adj, he, V);

    
    
    const double h = igl::avg_edge_length(V, F);
    double t = pow(h, 2);

    Eigen::MatrixXd U0=Eigen::MatrixXd::Zero(V.rows(), 1);
    Eigen::MatrixXd U1 = Eigen::MatrixXd::Zero(V.rows(), 1);

    //U0(122, 0) = 1;
    U0(10, 0) = 1;
    

    
    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    
    
    A = M - t * L;
    
    Eigen::FullPivLU<Eigen::MatrixXd> dec(A);
    U1 = dec.solve(U0);

    
    // Compute gradient operator:
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);

    
    
    
    Eigen::MatrixXd GU = Eigen::Map<const Eigen::MatrixXd>((G * U1).eval().data(), F.rows(), 3);
    
    MatrixXd div;
    MatrixXd area;

    area = areaMatrix(V, F);
    std::cout << G.rows() << " " << G.cols() << std::endl;


  
    
    //div=divMatrix1(adj, V,F, area,MatrixXd(M),MatrixXd(G));
    
    
    GU = -GU;
    GU.rowwise().normalize();
    div = divMatrix(adj, V, F, area, MatrixXd(M), GU);

    //flatettening the Gradient matrix
    Eigen::MatrixXd GUFlatttend;
    GUFlatttend = GU.transpose();  
    GUFlatttend.resize(GUFlatttend.rows() * GUFlatttend.cols(), 1);  
    
    MatrixXd phi=MatrixXd::Zero(V.rows(),1);

    Eigen::FullPivLU<Eigen::MatrixXd> dec1(L);
    //phi=dec1.solve(div*GUFlatttend);

    //phi = dec1.solve(div);
      
    //showVectorField(V, F, GU, viewer);
   
    
    // Initialize white
    C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
    //viewer.data().set_colors(C);
    

    viewer.callback_mouse_down = [&V, &F, &C,&he,&phi,&U1](igl::opengl::glfw::Viewer& viewer, int, int)->bool
    {
        int fid;
        Eigen::Vector3f bc;
        // Cast a ray in the view direction starting from the mouse position
        double x = viewer.current_mouse_x;
        
       
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;

        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
            viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
        {
            // paint hit red
            C.row(fid) << 1, 0, 0;
            
            std::cout << "index is " << F.row(fid) << std::endl;
            
            std::cout << "index is " << phi(F(fid,0),0) << std::endl;
            std::cout << "index is " << phi(F(fid, 1),0) << std::endl;
            std::cout << "index is " << phi(F(fid, 2), 0) << std::endl;


            for (const auto& value : vertexNeighbour(he, F(fid, 0))) {
                C.row(value) << 1, 0, 0;
            }
            for (const auto& value : vertexNeighbour(he, F(fid, 1))) {
                C.row(value) << 0, 1, 0;
            }
            for (const auto& value : vertexNeighbour(he, F(fid, 2))) {
                C.row(value) << 0, 0, 1;
            }
            //viewer.data().set_colors(C);
            return true;
        }
        return false;
    };
    
    
    
    std::cout << R"(Usage:[click]  Pick face on shape)" << std::endl;
    // Show mesh
    viewer.data().set_mesh(V, F);

    
    
    
    //std::cout << div << std::endl;
   
    
    for (int i = 0; i < phi.rows(); i = i + 1) {
         phi(i, 0) +=-phi.minCoeff() ;
    }

    for (int i = 0; i < phi.rows(); i = i + 1) {
        phi(i, 0) = floor(phi(i, 0));
            
        
    }

    showMaxMinCurvature(V, F, viewer);
    

    
    viewer.data().set_data(div);   

    viewer.callback_key_down = &key_down;
    viewer.data().show_lines = true;
    viewer.launch();
    
}










