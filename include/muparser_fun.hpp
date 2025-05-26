#ifndef _MUPARSERFUN_HPP_
#define _MUPARSERFUN_HPP_

#include "params.hpp"
#include <muParser.h>

#include <memory>
#include <string>

namespace laplaceSolver {
    class MuparserFun
    {
    public:
    MuparserFun(const MuparserFun &m)
        : m_parser(m.m_parser)
    {   try{
            m_parser.DefineVar("x1", &m_x1);
            m_parser.DefineVar("x2", &m_x2);
            m_parser.DefineConst("pi", pi);
        }
        catch(mu::Parser::exception_type &e){
            std::cerr << e.GetMsg() << std::endl;
        }

    };

    MuparserFun(): m_parser(){
         try{
            m_parser.DefineVar("x1", &m_x1);
            m_parser.DefineVar("x2", &m_x2);
            m_parser.DefineConst("pi", pi);
}
        catch(mu::Parser::exception_type &e){
            std::cerr << e.GetMsg() << std::endl;
        }

    }

    MuparserFun(const std::string &s)
    {
        try
        {
            m_parser.DefineVar("x1", &m_x1);
            m_parser.DefineVar("x2", &m_x2);
            m_parser.DefineConst("pi", pi);

            m_parser.SetExpr(s);
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cerr << e.GetMsg() << std::endl;
        }
    };

    scalar_type
    operator()(const node_type & x)
    {
        m_x1 = x[0];
        m_x2 = x[1];
        scalar_type y = m_parser.Eval();
        return y;
    };

    scalar_type operator()(const scalar_type x1, const scalar_type x2)
    {
        m_x1 = x1;
        m_x2 = x2;
        scalar_type y = m_parser.Eval();
        return y;
    };


    void SetExpression(const std::string & expr){
    try{
        m_parser.SetExpr(expr);
    }
    catch(mu::Parser::exception_type &e){
        std::cerr << e.GetMsg() << std::endl;
        }
    }

    private:
    scalar_type     m_x1;
    scalar_type     m_x2;

    mu::Parser m_parser;
    };

} // namespace laplaceSolver

#endif // _MUPARSERFUN_HPP_
