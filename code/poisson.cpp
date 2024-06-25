//
// Created by Liogky Alexey on 05.04.2022.
//

/**
 * Here we show the simplest way to write code used the P_1 templated space, 
 * the code essentially uses the structure of P_1 space and does not allow you to change the space without writing additional code 
 * 
 * This program generates and solves a finite element system for the stationary diffusion problem
 * \f[
 * \begin{aligned}
 *  -\mathrm{div}\ \mathbf{D}\ &\mathrm{grad}\ u\ &= \mathbf{F}\ in\  \Omega          \\
 *                             &               u\ &= U_0\        on\  \partial \Omega \\
 * \end{aligned}
 * \f]
 * where Ω = [0,1]^3
 * The user-defined coefficients are
 *    D(x)   = I           - identity tensor
 *    F(x)   = 1           - right-hand side
 *    U_0(x) = 0           - essential (Dirichlet) boundary condition
 *
 * @see src/Tutorials/PackageFEM/main_simple.f in [Ani3d library](https://sourceforge.net/projects/ani3d/)
 */
#include "inmost.h"
#include "mesh_gen.h"
#include "anifem++/inmost_interface/fem.h"

using namespace INMOST;

int main(){

    // Инициализация сетки и решателя
    int pRank = 0, pCount = 1;
    Mesh::Initialize(0, NULL);
    Solver::Initialize(0, NULL);

    // Генерация структурированной равномерной сетки куба 
    std::unique_ptr<Mesh> mptr = GenerateCube(INMOST_MPI_COMM_WORLD, 8, 8, 8);
    Mesh* m = mptr.get();
    // Вывод параметров сетки
    long nN = m->TotalNumberOf(NODE), nE = m->TotalNumberOf(EDGE), nF = m->TotalNumberOf(FACE), nC = m->TotalNumberOf(CELL);
    std::cout << "Mesh info:" << " #N " << nN << " #E " << nE << " #F " << nF << " #T " << nC << std::endl;
    using namespace Ani;
    // Выбор пространства конечных элементов
    using UFem = FemFix<FEM_P1>;
    constexpr auto UNF = Operator<IDEN, UFem>::Nfa::value;
    // Определение меток границ
    Tag BndLabel = m->CreateTag("bnd_label", DATA_INTEGER, NODE, NONE, 1);
    for (auto it = m->BeginNode(); it != m->EndNode(); ++it){
        int lbl = 0;
        for (auto k = 0; k < 3; ++k)
            if (abs(it->Coords()[k] - 0) < 10*std::numeric_limits<double>::epsilon())
                lbl |= (1 << (2*k+0));
            else if (abs(it->Coords()[k] - 1) < 10*std::numeric_limits<double>::epsilon()) 
                lbl |= (1 << (2*k+1));   
        it->Integer(BndLabel) = lbl;
    }

    // Определение тензоров задачи
    auto D_tensor = [](const Coord<> &X, double *D, TensorDims dims, void *user_data, int iTet) {
        (void) dims; (void) user_data; (void) iTet;
        D[0] = 1;  
        return Ani::TENSOR_SCALAR;
    };
    auto F_tensor = [](const Coord<> &X, double *F, TensorDims dims, void *user_data, int iTet) {
        (void) X; (void) dims; (void) user_data; (void) iTet;
        F[0] = 1;
        return Ani::TENSOR_SCALAR;
    };
    auto U_0 = [](const Coord<> &X, double* res, ulong dim, void* user_data)->int{ (void) X; (void) dim; (void) user_data; res[0] = 0.0; return 0; }; 

    // Определение тега для записи решения
    Tag u = m->CreateTag("u", DATA_REAL, NODE, NONE, 1);

    // Определение пользовательской структуры для ввода пользовательских данных в локальный сборщик
    struct ProbLocData{
        // Сохранение меток используемых для задания граничных условий 
        std::array<int, 4> nlbl = {0};
    };
 
    // Определение элементарного сборщика локальных матриц и правых частей
    std::function<void(const double**, double*, double*, void*)> local_assembler =
            [&D_tensor, &F_tensor, &U_0, order = 2](const double** XY/*[4]*/, double* Adat, double* Fdat, void* user_data) -> void{
        DenseMatrix<> A(Adat, UNF, UNF), F(Fdat, UNF, 1);
        A.SetZero(); F.SetZero();

        Tetra<const double> XYZ(XY[0], XY[1], XY[2], XY[3]);
        // Здесь мы задаем DfuncTraits<TENSOR_SCALAR, true> потому что знаем, что тензор является скаляром и постоянной
        // Элементарная матрица жесткости <grad(P1), D grad(P1)>
        fem3Dtet<Operator<GRAD, UFem>, Operator<GRAD, UFem>, DfuncTraits<TENSOR_SCALAR, true>>(XYZ, D_tensor, A, order);

        // Элементарный вектор правой части <F, P1>
        fem3Dtet<Operator<IDEN, FemFix<FEM_P0>>, Operator<IDEN, UFem>, DfuncTraits<TENSOR_SCALAR, true>>(XYZ, F_tensor, F, order);

        // Считываем метки вершин из пользовательских данных
        auto& dat = *static_cast<ProbLocData*>(user_data);
        // Применяем граничные условия Дирихле
        for (int i = 0; i < 4; ++i)
            if (dat.nlbl[i] > 0) {
                double bc = 0;
                U_0(std::array<double, 3>{XY[i][0], XY[i][1], XY[i][2]}, &bc, 1, nullptr);
                applyDir(A, F, i, bc); // задание граничного условия Дирихле
            }    
    };
    // Определение функции собирающей данные из каждого тетраэдра для передачи из в локальный сборщик
    auto local_data_gatherer = [&BndLabel](ElementalAssembler& p) -> void{
        double *nn_p = p.get_nodes();
        const double *args[] = {nn_p, nn_p + 3, nn_p + 6, nn_p + 9};
        ProbLocData data;
        std::fill(data.nlbl.begin(), data.nlbl.end(), 0);
        for (unsigned i = 0; i < data.nlbl.size(); ++i) 
            data.nlbl[i] = (*p.nodes)[i].Integer(BndLabel);    

        p.compute(args, &data);
    };

    // Определение сборщика
    Assembler discr(m);
    discr.SetMatRHSFunc(GenerateElemMatRhs(local_assembler, UNF, UNF));
    {
        // Создание счетчика глобальных степеней свободы
        auto Var0Helper = GenerateHelper<UFem>(); //DofT::DofMap(Dof<UFem>::Map());
        FemExprDescr fed;
        fed.PushTrialFunc(Var0Helper, "u");
        fed.PushTestFunc(Var0Helper, "phi_u");
        discr.SetProbDescr(std::move(fed));
    }
    discr.SetDataGatherer(local_data_gatherer);
    discr.PrepareProblem();
    Sparse::Matrix A("A");
    Sparse::Vector x("x"), b("b");
    // Инициализация матрицы и правой части (является опциональной)
    A.SetInterval(discr.getBegInd(), discr.getEndInd()); // Инициализация в виде свободной матрицы
    b.SetInterval(discr.getBegInd(), discr.getEndInd()); // Инициализация нулевым вектором
    // В этой задаче нам не требуется информация по граням и ребрам, поэтому мы отключаем сбор данных по ним
    discr.m_assm_traits.prepare_faces = discr.m_assm_traits.prepare_edges = false;
    // Сборка матрицы и правой части
    TimerWrap m_timer_total; m_timer_total.reset();
    discr.Assemble(A, b);
    double total_assm_time = m_timer_total.elapsed();
    // Вывод времени сборки
    std::cout << "#dofs = " << discr.m_enum.getMatrixSize() << std::endl; 
    std::cout << "Assembler timers:" << "\n\tTotal         : " << total_assm_time << "s" << std::endl;

    // Настройка линейного решателя и решение собранной системы
    Solver solver("inner_mptiluc");
    solver.SetMatrix(A);
    solver.Solve(b, x);
    
    // Копирование результатов в тэг и сохранение решения
    discr.SaveVar(x, 0, u);
    m->Save("poisson.pvtu");

    discr.Clear();
    mptr.reset();
    Solver::Finalize();
    Mesh::Finalize();

    return 0;
}
