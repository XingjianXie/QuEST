// read xml file from all.xml in js
import fs from 'fs';
import util from 'util';
import { XMLParser } from 'fast-xml-parser';

function defmapper(def) {
    return def.map(member => {
        if (member.references && !Array.isArray(member.references)) {
            member.references = [member.references];
        }
        return {
            name: member.name,
            references: member.references ? member.references.filter(
                ref =>
                ref['@_refid'].startsWith('QuEST__common_8c') ||
                ref['@_refid'].startsWith('QuEST__internal_8h')
                ).map(
                    ref => {
                        return {
                            name: ref['#text'],
                            type: ref['@_refid'].startsWith('QuEST__common_8c') ? 'common' : 'internal'
                        }
                    }
                ) : []
        }
    })
}

function addref(common, gf, fr, set) {
    for (const ref of fr) {
        if (ref.type === 'internal') {
            gf.push({name: ref.name, type: 'internal'});
        } else {
            const elem = {name: ref.name, type: 'common', dependencies: []};
            if (set.has(ref.name)) {
                elem.dependencies = 'recursive'
                gf.push(elem);
            } else {
                set.add(ref.name);
                addref(common, elem.dependencies, common.find(c => c.name === ref.name)?.references ?? [], set)
                gf.push(elem);
                set.delete(ref.name);
            }
        }
    }
}

async function main() {
    const xml = fs.readFileSync('./doc/xml/all.xml', 'utf-8');
    const parser = new XMLParser({
        ignoreAttributes: false,
    });
    const obj = parser.parse(xml);
    const compounddef = obj.doxygen.compounddef.reduce((acc, cur) => {
        acc[cur.compoundname] = cur
        return acc
    }, {});
    const unitarydef = compounddef.unitary.sectiondef.memberdef;
    const unitary = defmapper(unitarydef);
    const normgatedef = compounddef.normgate.sectiondef.memberdef;
    const normgate = defmapper(normgatedef);
    const operatordef = compounddef.operator.sectiondef.memberdef;
    const operator = defmapper(operatordef);
    const commondef = compounddef['QuEST_common.c'].sectiondef.find(def => def['@_kind'] === 'func').memberdef;
    const common = defmapper(commondef);

    const graph = [];

    for (const func of unitary) {
        const elem = {name: func.name, type: 'unitary', dependencies: []};
        addref(common, elem.dependencies, func.references, new Set());
        graph.push(elem);
    }

    for (const func of normgate) {
        const elem = {name: func.name, type: 'normgate', dependencies: []};
        addref(common, elem.dependencies, func.references, new Set());
        graph.push(elem);
    }

    for (const func of operator) {
        const elem = {name: func.name, type: 'operator', dependencies: []};
        addref(common, elem.dependencies, func.references, new Set());
        graph.push(elem);
    }

    console.log(JSON.stringify(graph, null, 2));
}

main();
