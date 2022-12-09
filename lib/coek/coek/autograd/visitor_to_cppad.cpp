#include <unordered_map>

#include "../ast/base_terms.hpp"
#include "../ast/constraint_terms.hpp"
#include "../ast/expr_terms.hpp"
#include "../ast/value_terms.hpp"
#include "coek/api/constraint.hpp"
#include "coek/api/objective.hpp"
#include "coek/model/model.hpp"
#include "coek/model/model_repn.hpp"
#include "cppad_repn.hpp"

namespace coek {

//
// This empty namespace contains functions used to walk the COEK
// expression tree.  The CppAD expression is accumulated in the 'ans'
// value.
//
namespace {

class VisitorData {
   public:
    std::unordered_map<expr_pointer_t, CppAD::AD<double> > cache;

    std::vector<CppAD::AD<double> >& ADvars;
    std::unordered_map<VariableRepn, size_t>& used_variables;
    std::map<VariableRepn, size_t>& fixed_variables;
    std::map<ParameterRepn, size_t>& parameters;
    std::vector<CppAD::AD<double> >& dynamic_params;

    VisitorData(std::vector<CppAD::AD<double> >& _ADvars,
                std::unordered_map<VariableRepn, size_t>& _used_variables,
                std::map<VariableRepn, size_t>& _fixed_variables,
                std::map<ParameterRepn, size_t>& _parameters,
                std::vector<CppAD::AD<double> >& _dynamic_params)
        : ADvars(_ADvars),
          used_variables(_used_variables),
          fixed_variables(_fixed_variables),
          parameters(_parameters),
          dynamic_params(_dynamic_params)
    {
    }
};

void visit_expression(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans);

void visit_ConstantTerm(const expr_pointer_t& expr, VisitorData& /*data*/, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<ConstantTerm>(expr);
    ans += tmp->value;
}

void visit_ParameterTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<ParameterTerm>(expr);
    ans += data.dynamic_params[data.parameters[tmp]];
}

void visit_VariableTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<VariableTerm>(expr);
    if (tmp->fixed)
        ans += data.dynamic_params[data.fixed_variables[tmp]];
    else
        ans += data.ADvars[data.used_variables[tmp]];
}

void visit_MonomialTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<MonomialTerm>(expr);
    if (tmp->var->fixed)
        ans += tmp->coef * data.dynamic_params[data.fixed_variables[tmp->var]];
    else
        ans += tmp->coef * data.ADvars[data.used_variables[tmp->var]];
}

#define FROM_BODY(TERM)                                                                \
    void visit_##TERM(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans) \
    {                                                                                  \
        auto tmp = std::dynamic_pointer_cast<TERM>(expr);                              \
        visit_expression(tmp->body, data, ans);                                        \
    }

FROM_BODY(InequalityTerm)
FROM_BODY(EqualityTerm)
FROM_BODY(ObjectiveTerm)

void visit_NegateTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<NegateTerm>(expr);
    CppAD::AD<double> body;
    visit_expression(tmp->body, data, body);
    ans += -body;
}

void visit_PlusTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<PlusTerm>(expr);

    for (auto& it : *(tmp->data)) {
        CppAD::AD<double> next;
        visit_expression(it, data, next);
        ans += next;
    }
}

void visit_TimesTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<TimesTerm>(expr);
    CppAD::AD<double> lhs;
    visit_expression(tmp->lhs, data, lhs);
    CppAD::AD<double> rhs;
    visit_expression(tmp->rhs, data, rhs);
    ans += lhs * rhs;
}

void visit_DivideTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<DivideTerm>(expr);
    CppAD::AD<double> lhs;
    visit_expression(tmp->lhs, data, lhs);
    CppAD::AD<double> rhs;
    visit_expression(tmp->rhs, data, rhs);
    ans += lhs / rhs;
}

#define UNARY_VISITOR(TERM, FN)                                                        \
    void visit_##TERM(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans) \
    {                                                                                  \
        auto tmp = std::dynamic_pointer_cast<TERM>(expr);                              \
        CppAD::AD<double> body;                                                        \
        visit_expression(tmp->body, data, body);                                       \
        ans += CppAD::FN(body);                                                        \
    }

UNARY_VISITOR(AbsTerm, abs)
// UNARY_VISITOR(CeilTerm, ceil)
// UNARY_VISITOR(FloorTerm, floor)
UNARY_VISITOR(ExpTerm, exp)
UNARY_VISITOR(LogTerm, log)
UNARY_VISITOR(Log10Term, log10)
UNARY_VISITOR(SqrtTerm, sqrt)
UNARY_VISITOR(SinTerm, sin)
UNARY_VISITOR(CosTerm, cos)
UNARY_VISITOR(TanTerm, tan)
UNARY_VISITOR(SinhTerm, sinh)
UNARY_VISITOR(CoshTerm, cosh)
UNARY_VISITOR(TanhTerm, tanh)
UNARY_VISITOR(ASinTerm, asin)
UNARY_VISITOR(ACosTerm, acos)
UNARY_VISITOR(ATanTerm, atan)
UNARY_VISITOR(ASinhTerm, asinh)
UNARY_VISITOR(ACoshTerm, acosh)
UNARY_VISITOR(ATanhTerm, atanh)

void visit_PowTerm(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto tmp = std::dynamic_pointer_cast<PowTerm>(expr);
    CppAD::AD<double> lhs;
    visit_expression(tmp->lhs, data, lhs);
    if (tmp->rhs->is_constant()) {
        double val = tmp->rhs->eval();
        if (fabs(val - int(val)) < 1e-12) {
            ans += CppAD::pow(lhs, int(val));
            return;
        }
    }

    CppAD::AD<double> rhs;
    visit_expression(tmp->rhs, data, rhs);
    ans += CppAD::pow(lhs, rhs);
}

#define VISIT_CASE(TERM)               \
    case TERM##_id: {                  \
        visit_##TERM(expr, data, ans); \
        data.cache[expr] = ans;        \
    } break

void visit_expression(const expr_pointer_t& expr, VisitorData& data, CppAD::AD<double>& ans)
{
    auto curr = data.cache.find(expr);
    if (curr != data.cache.end()) {
        ans += curr->second;
        return;
    }

    switch (expr->id()) {
        VISIT_CASE(ConstantTerm);
        VISIT_CASE(ParameterTerm);
        VISIT_CASE(VariableTerm);
        VISIT_CASE(MonomialTerm);
        VISIT_CASE(InequalityTerm);
        VISIT_CASE(EqualityTerm);
        VISIT_CASE(ObjectiveTerm);
        VISIT_CASE(NegateTerm);
        VISIT_CASE(PlusTerm);
        VISIT_CASE(TimesTerm);
        VISIT_CASE(DivideTerm);
        VISIT_CASE(AbsTerm);
        // VISIT_CASE(CeilTerm);
        // VISIT_CASE(FloorTerm);
        VISIT_CASE(ExpTerm);
        VISIT_CASE(LogTerm);
        VISIT_CASE(Log10Term);
        VISIT_CASE(SqrtTerm);
        VISIT_CASE(SinTerm);
        VISIT_CASE(CosTerm);
        VISIT_CASE(TanTerm);
        VISIT_CASE(SinhTerm);
        VISIT_CASE(CoshTerm);
        VISIT_CASE(TanhTerm);
        VISIT_CASE(ASinTerm);
        VISIT_CASE(ACosTerm);
        VISIT_CASE(ATanTerm);
        VISIT_CASE(ASinhTerm);
        VISIT_CASE(ACoshTerm);
        VISIT_CASE(ATanhTerm);
        VISIT_CASE(PowTerm);

        default:
            throw std::runtime_error(
                "Error in CppAD_Repn visitor!  Visiting unexpected expression term "
                + std::to_string(expr->id()));
    };
}

}  // namespace

void build_expression(const expr_pointer_t& root, std::vector<CppAD::AD<double> >& ADvars,
                                  CppAD::AD<double>& ans,
                                  std::unordered_map<VariableRepn, size_t>& _used_variables,
                                  std::map<VariableRepn, size_t>& fixed_variables,
                                  std::map<ParameterRepn, size_t>& parameters,
                                  std::vector<CppAD::AD<double> >& dynamic_params)
{
    VisitorData data(ADvars, _used_variables, fixed_variables, parameters, dynamic_params);
    visit_expression(root, data, ans);
}

}  // namespace coek
